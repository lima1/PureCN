globalVariables(names=c("Gene", "LR", "chrom", "seg.id", "seg.length",
"seg.mean"))

.calcLlikSegment <- function(subclonal, lr, sd.seg, p, Ci, total.ploidy, 
max.exon.ratio) {
    if (subclonal) {
        return(.calcLlikSegmentSubClonal(lr, max.exon.ratio))
    }

    .calcLlikSegmentClonal(lr, sd.seg, p, Ci, total.ploidy)
}
.calcLlikSegmentClonal <- function(lr, sd.seg, p, Ci, total.ploidy) {
    sum(dnorm(lr, mean = log2((p * Ci + (1 - p) * 2)/total.ploidy), 
        sd = sd.seg, log = TRUE))
}
.calcLlikSegmentSubClonal <- function(lr, max.exon.ratio) {
    sum(dunif(x = sapply(2^lr, function(y) min(y, max.exon.ratio)), min = 0, 
        max = max.exon.ratio, log = TRUE))
}
.calcLogRatio <- function(normal, tumor, verbose) {
    # make sure that normal and tumor align
    if (!identical(as.character(normal[, 1]), as.character(tumor[, 1]))) {
        stop("Interval files in normal and tumor different.")
    }
    total.cov.normal <- sum(as.numeric(normal$coverage), na.rm = TRUE)
    total.cov.tumor <- sum(as.numeric(tumor$coverage), na.rm = TRUE)

    log.ratio <- log2(tumor$average.coverage/normal$average.coverage) + 
                 log2(total.cov.normal/total.cov.tumor)

    mean.log.ratio <- mean(subset(log.ratio, !is.infinite(log.ratio)), 
        na.rm = TRUE)
    # calibrate
    log.ratio <- log.ratio - mean.log.ratio
    log.ratio
}
.calcSNVLLik <- function(vcf, tumor.id.in.vcf, ov, p, test.num.copy, 
    C.posterior, C, snv.model, prior.somatic, snv.lr, sampleid = NULL, 
    cont.rate = 0.01, prior.M = NULL, post.optimize) {

    prior.cont <- ifelse(info(vcf)$DB, cont.rate, 0)
    prior.somatic <- prior.somatic/(1 + cont.rate)
    
    haploid.penalty <- 0
    
    if (median(C) < 1.1 && p < 0.3) {
        haploid.penalty <- 1
    }
    
    subclonal <- apply(C.posterior[queryHits(ov), ], 1, which.max) == ncol(C.posterior)
    
    # very high amplifications results in prior probabilities for tested copy numbers
    # of 0. We add the subclonal posterior to the maximum copy number to fix this.
    idx <- abs(round(C) - C) > 0.001 & C > max(test.num.copy)
    if (length(idx) > 1) {
        C.posterior[idx, ncol(C.posterior) - 1] <- C.posterior[idx, 
            ncol(C.posterior) - 1] + C.posterior[idx, ncol(C.posterior)]
    }
    
    seg.idx <- which(1:nrow(C.posterior) %in% queryHits(ov))
    sd.ar <- sd(unlist(geno(vcf)$FA[, tumor.id.in.vcf]))
    xx <- lapply(seg.idx, function(i) {
        # classify germline vs somatic
        idx <- subjectHits(ov)[queryHits(ov) == i]
        
        ar_all <- unlist(geno(vcf)$FA[idx, tumor.id.in.vcf])
        ar_all[ar_all > 1] <- 1
        dp_all <- unlist(geno(vcf)$DP[idx, tumor.id.in.vcf])

        lapply(test.num.copy, function(Ci) {
            lr.C <- rep(Ci, length(ar_all))
            pM <- prior.M[[as.character(i)]]
            pM.both <- list(rep(log(Ci + 1 + haploid.penalty), 
                                length(test.num.copy)), 
                            rep(log(Ci + 1 + haploid.penalty), 
                                length(test.num.copy)))
            
            if (!is.null(pM)) {
                pM.both <- list(rep(log(Ci + 1 + haploid.penalty), 
                                    length(test.num.copy)), 
                                log((1 - pM) + (Ci + 1) + haploid.penalty))
            }
            
            p.ar <- lapply(c(0, 1), function(g) 
                lapply(1:length(ar_all), function(j) 
                sapply(test.num.copy, function(Mi) 
                    ifelse(Mi > lr.C[j], -Inf, 
                    dbeta(x = (p * Mi + g * (1-p))/(p * lr.C[j] + 2 * (1-p)), 
                    shape1 = ar_all[j] * dp_all[j] + 1, 
                    shape2 = (1 - ar_all[j]) * dp_all[j] + 1, log = TRUE) -
                    pM.both[[g + 1]][which.min(abs(Mi - test.num.copy))]))))
            
            p.ar.cont.1 <- lapply(1:length(ar_all), function(j) 
                dbeta(x = (p * lr.C[j] + 2 * (1 - p - cont.rate))/
                    (p * lr.C[j] + 2 * (1 - p)), 
                    shape1 = ar_all[j] * dp_all[j] + 1, 
                    shape2 = (1 - ar_all[j]) * dp_all[j] + 1, log = TRUE) - 
                log(lr.C[j] + 1 + haploid.penalty))

            p.ar.cont.2 <- lapply(1:length(ar_all), function(j) 
                dbeta(x = (cont.rate)/(p * lr.C[j] + 2 * (1 - p)), 
                    shape1 = ar_all[j] * dp_all[j] + 1, 
                    shape2 = (1 - ar_all[j]) * dp_all[j] + 1, log = TRUE) - 
                log(lr.C[j] + 1 + haploid.penalty))
            
            # add prior probabilities for somatic vs germline
            p.ar[[1]] <- lapply(1:length(p.ar[[1]]), 
                function(j) p.ar[[1]][[j]] + log(prior.somatic[idx][j]))

            p.ar[[2]] <- lapply(1:length(p.ar[[2]]), 
                function(j) p.ar[[2]][[j]] + log(1 - prior.somatic[idx][j]))

            # contamination (either homozygous germline, or germline from 
            # other sample)

            p.ar[[3]] <- lapply(1:length(p.ar.cont.1), 
                function(j) p.ar.cont.1[[j]] + log(prior.cont[idx][j]))
            p.ar[[4]] <- lapply(1:length(p.ar.cont.2), 
                function(j) p.ar.cont.2[[j]] + log(prior.cont[idx][j]))
            
            do.call(cbind, lapply(p.ar, function(x) do.call(rbind, x)))
        })
        
    })
    snv.posteriors <- do.call(rbind, 
        lapply(1:length(xx), function(i) Reduce("+", 
            lapply(test.num.copy, function(Ci) 
                exp(xx[[i]][[Ci + 1]]) * C.posterior[seg.idx[i], Ci + 1]))))

    colnames(snv.posteriors) <- c(paste("SOMATIC.M", test.num.copy, sep = ""), 
        paste("GERMLINE.M", test.num.copy, sep = ""), "GERMLINE.CONTHIGH", 
        "GERMLINE.CONTLOW")
    
    vcf.ids <- do.call(c, lapply(seg.idx, function(i) 
        subjectHits(ov)[queryHits(ov) == i]))
    rownames(snv.posteriors) <- vcf.ids
    
    # this just adds a lot of helpful info to the SNV posteriors
    xx <- .extractMLSNVState(snv.posteriors)
    
    posteriors <- snv.posteriors/rowSums(snv.posteriors)
    posteriors <- cbind(
        as.data.frame(rowRanges(vcf[vcf.ids]))[, 1:3], 
        posteriors, 
        xx, 
        ML.C = C[queryHits(ov)]
    )
    
    posteriors$ML.AR <- (p * posteriors$ML.M + 
        ifelse(posteriors$ML.SOMATIC, 0, 1) * 
        (1 - p))/(p * posteriors$ML.C + 2 * (1 - p))

    posteriors$AR <- unlist(geno(vcf[vcf.ids])$FA[, tumor.id.in.vcf])
    posteriors$CN.Subclonal <- subclonal
    posteriors$Log.Ratio <- snv.lr[vcf.ids]
    posteriors$Prior.Somatic <- prior.somatic[vcf.ids]
    posteriors$Prior.Contamination <- prior.cont[vcf.ids]
    
    # Extract LOH
    posteriors$ML.LOH <- (posteriors$ML.M == posteriors$ML.C | 
        posteriors$ML.M == 0 | posteriors$ML.C == 1)

    loh <- segment(
        CNA(ifelse(posteriors$ML.LOH, 1, 0), 
            chrom = posteriors$seqnames, 
            maploc = posteriors$start, 
            data.type = "binary", 
            sampleid = sampleid, 
            presorted = TRUE))
    
    # these are potential artifacts with very high clonal probability and would have
    # huge impact on log-likelihood
    rm.snv.posteriors <- apply(snv.posteriors, 1, max)
    idx.ignore <- (posteriors$CN.Subclonal & 
        posteriors$ML.C < max(test.num.copy) & 
        C.posterior[queryHits(ov), ncol(C.posterior)] > 0.95) | 
        rm.snv.posteriors == 0

    ret <- list(
        llik = sum(log(rm.snv.posteriors[!idx.ignore])) - sum(idx.ignore), 
        likelihoods = snv.posteriors, 
        posteriors = posteriors, 
        vcf.ids = vcf.ids, 
        segment.ids = queryHits(ov), 
        loh = loh, 
        llik.ignored = idx.ignore)

    if (post.optimize && is.null(prior.M)) {
        ret <- .calcSNVLLik(vcf, tumor.id.in.vcf, ov, p, test.num.copy, 
            C.posterior, C, snv.model, prior.somatic, snv.lr, sampleid, 
            cont.rate, prior.M = .calcMpriorGermline(ret), 
            post.optimize = post.optimize)
    }
    ret
}

.extractMLSNVState <- function(snv.posteriors) {
    l1 <- apply(snv.posteriors, 1, which.max)
    xx <- do.call(rbind, strsplit(colnames(snv.posteriors)[l1], "\\."))
    xx[, 1] <- ifelse(xx[, 1] == "GERMLINE", FALSE, TRUE)
    xx[, 2] <- gsub("^M", "", xx[, 2])
    xx <- as.data.frame(xx, stringsAsFactors = FALSE)
    xx[, 1] <- as.logical(xx[, 1])
    xx[, 2] <- suppressWarnings(as.numeric(xx[, 2]))
    colnames(xx) <- c("ML.SOMATIC", "ML.M")
    xx
}

.checkParameters <- function(test.purity, min.ploidy, max.ploidy, max.non.clonal) {
    if (min(test.purity) <= 0 || max(test.purity) > 1) 
        stop("test.purity not within expected range.")
    if (min.ploidy <= 0 || max.ploidy <= 2) 
        stop("min.ploidy or max.ploidy not within expected range.")
    if (max.non.clonal > 1) 
        stop("max.non.clonal not within expected range.")
}
.failedNonAberrant <- function(result, cutoffs = c(0.005, 0.002)) {
    xx <- split(result$seg, result$seg$C)
    if (length(xx) < 3) 
        return(TRUE)
    xx.sum <- sort(sapply(xx, function(x) sum(x$size)), decreasing = TRUE)
    xx.sum <- xx.sum/sum(xx.sum)
    if (xx.sum[2] <= cutoffs[1] && xx.sum[3] <= cutoffs[2]) 
        return(TRUE)
    FALSE
}
.filterDuplicatedResults <- function(results) {
    if (length(results) < 2) 
        return(results)
    idx <- 2:length(results)
    diff.purity <- abs(sapply(results[idx - 1], function(x) x$purity) - sapply(results[idx], 
        function(x) x$purity))
    diff.ploidy <- abs(sapply(results[idx - 1], function(x) x$ploidy) - sapply(results[idx], 
        function(x) x$ploidy))
    idx.duplicated <- c(FALSE, diff.purity < 0.1 & diff.ploidy < 0.2)
    results[!idx.duplicated]
}
.findLocalMinima <- function(m) {
    loc.min <- matrix(nrow = 0, ncol = 2)
    for (i in 1:nrow(m)) {
        for (j in 1:ncol(m)) {
            x <- (i - 1):(i + 1)
            x <- x[x >= 1 & x <= nrow(m)]
            y <- (j - 1):(j + 1)
            y <- y[y >= 1 & y <= ncol(m)]
            if (m[i, j] == max(m[x, y]) && !is.infinite(m[i, j])) 
                loc.min <- rbind(loc.min, c(row = i, col = j))
        }
    }
    loc.min
}
.appendComment <- function(a, b) {
    if (is.na(a)) 
        return(b)
    paste(a, b, sep = ";")
}
.flagResult <- function(result, max.non.clonal = 0.2) {
    result$flag_comment <- NA
    result$flag <- .failedNonAberrant(result)
    if (result$flag) {
        result$flag_comment <- .appendComment(result$flag_comment, "NON-ABERRANT")
    }
    if (result$fraction.subclonal > max.non.clonal) {
        result$flag <- TRUE
        result$flag_comment <- .appendComment(result$flag_comment, "POLYGENOMIC")
    }
    if (result$ploidy > 4.5 || result$ploidy < 1.5) {
        result$flag <- TRUE
        result$flag_comment <- .appendComment(result$flag_comment, "RARE KARYOTYPE")
    }
    if (!is.null(result$SNV.posterior)) {
        tmp <- result$SNV.posterior$beta.model$loh$output
        if (!is.null(tmp)) {
            tmp <- tmp[complete.cases(tmp), ]
            fraction.loh <- sum(tmp[which(tmp$seg.mean > 0.2), "num.mark"])/sum(tmp$num.mark)
            # Assume that everything below 2.6 did not undergo genome duplication, which can
            # result in lots of LOH
            if (result$ploidy < 2.6 && fraction.loh > 0.5) {
                result$flag <- TRUE
                result$flag_comment <- .appendComment(result$flag_comment, "EXCESSIVE LOH")
            }
        }
    }
    return(result)
}

.flagResults <- function(results, max.non.clonal = 0.2, max.logr.sdev, logr.sdev, max.segments,
    flag = NA, flag_comment = NA) {
    results <- lapply(results, .flagResult, max.non.clonal = max.non.clonal)
    
    # ldiff <- ( results[[1]]$total.log.likelihood -
    # results[[2]]$total.log.likelihood ) / abs( results[[1]]$total.log.likelihood )
    # if (ldiff < 0.1) { results[[1]]$flag <- TRUE results[[1]]$flag_comment <-
    # 'SIMILAR TO SECOND MOST LIKELY OPTIMUM' return(results) }
    number.segments <- nrow(results[[1]]$seg)
    
    if (logr.sdev > max.logr.sdev) {
        for (i in 1:length(results)) {
            results[[i]]$flag <- TRUE
            results[[i]]$flag_comment <- .appendComment(results[[i]]$flag_comment, 
                "NOISY LOG-RATIO")
        }
    }

    if (number.segments > max.segments) {
        for (i in 1:length(results)) {
            results[[i]]$flag <- TRUE
            results[[i]]$flag_comment <- .appendComment(results[[i]]$flag_comment, 
                "NOISY SEGMENTATION")
        }
    }
    
    # some global flags
    if (!is.na(flag) && flag) {
        for (i in 1:length(results)) {
            results[[i]]$flag <- TRUE
            results[[i]]$flag_comment <- .appendComment(results[[i]]$flag_comment, 
                flag_comment)
        }
    }
    results
}
.getGeneCalls <- function(seg.adjusted, gc.data, log.ratio, fun.focal, args.focal) {
    args.focal <- c(list(seg = seg.adjusted), args.focal)
    focal <- do.call(fun.focal, args.focal)
    abs.gc <- GRanges(seqnames = seg.adjusted$chrom, IRanges(start = seg.adjusted$loc.start, 
        end = seg.adjusted$loc.end))
    
    gc.pos <- as.data.frame(do.call(rbind, strsplit(gc.data$Target, ":|-")), stringsAsFactors = FALSE)
    colnames(gc.pos) <- c("chrom", "start", "end")
    gc.pos <- cbind(Gene = gc.data$Gene, gc.pos)
    gc.pos$start <- as.numeric(gc.pos$start)
    gc.pos$end <- as.numeric(gc.pos$end)
    
    # that will otherwise mess up the log-ratio means, mins and maxs
    idx <- !is.nan(log.ratio) & !is.infinite(log.ratio)
    gc.pos <- gc.pos[idx, ]
    log.ratio <- log.ratio[idx]
    
    gc.gr <- GRanges(seqnames = .strip.chr.name(gc.pos$chrom), IRanges(start = gc.pos$start, 
        end = gc.pos$end))
    ov <- findOverlaps(gc.gr, abs.gc)
    # use funky data.table to calculate means etc. in two lines of code.
    dt <- data.table(gc.pos[queryHits(ov), ], C = seg.adjusted[subjectHits(ov), "C"], 
        seg.mean = seg.adjusted[subjectHits(ov), "seg.mean"], LR = log.ratio[queryHits(ov)], 
        seg.id = subjectHits(ov), seg.length = seg.adjusted$size[subjectHits(ov)], 
        focal = focal[subjectHits(ov)])
    gene.calls <- data.frame(dt[, list(chr = chrom[1], start = min(start), end = max(end), 
        C = median(C), seg.mean = median(seg.mean), seg.id = seg.id[which.min(seg.length)], 
        number.exons = length(start), gene.mean = mean(LR, na.rm = TRUE), gene.min = min(LR, 
            na.rm = TRUE), gene.max = max(LR, na.rm = TRUE), focal = focal[which.min(seg.length)]), 
        by = Gene], row.names = 1)
}

.optimizeGrid <- function(test.purity, min.ploidy, max.ploidy, test.num.copy = 0:7, 
    exon.lrs, seg, sd.seg, li, max.exon.ratio, max.non.clonal, verbose = FALSE, debug = FALSE) {
    ploidy.grid <- seq(min.ploidy, max.ploidy, by = 0.2)
    if (min.ploidy < 1.8 && max.ploidy > 2.2) {
        ploidy.grid <- c(seq(min.ploidy, 1.8, by = 0.2), 1.9, 2, 2.1, seq(2.2, max.ploidy, 
            by = 0.2))
    }
    mm <- lapply(test.purity, function(p) {
        b <- 2 * (1 - p)
        log.ratio.offset <- 0
        if (debug) 
            message(paste(b, log.ratio.offset))
        lapply(ploidy.grid, function(D) {
            dt <- p/D
            llik.all <- lapply(1:length(exon.lrs), function(i) .calcLlikSegmentExonLrs(exon.lrs[[i]], 
                log.ratio.offset, max.exon.ratio, sd.seg, dt, b, D, test.num.copy))
            subclonal <- sapply(llik.all, which.max) == 1
            subclonal.f <- length(unlist(exon.lrs[subclonal]))/length(unlist(exon.lrs))
            if (debug) 
                message(paste(sum(subclonal), subclonal.f))
            if (subclonal.f > max.non.clonal) 
                return(-Inf)
            llik <- sum(sapply(llik.all, max))
        })
    })
    mm <- sapply(mm, function(x) unlist(x))
    colnames(mm) <- test.purity
    rownames(mm) <- ploidy.grid
    ai <- .findLocalMinima(mm)
    candidates <- data.frame(ploidy = as.numeric(rownames(mm)[ai[, 1]]), purity = as.numeric(colnames(mm)[ai[, 
        2]]), llik = mm[ai])
    candidates$tumor.ploidy <- (candidates$ploidy - 2 * (1 - candidates$purity))/candidates$purity
    
    # add a diploid candidate solution in case there is none.
    if (min(abs(2 - candidates$tumor.ploidy)) > 0.3) {
        candidates <- rbind(candidates, c(2, as.numeric(names(which.max(mm["2", ]))), 
            max(mm["2", ]), 2))
    }
    
    candidates <- candidates[candidates$tumor.ploidy >= 0.5, ]
    
    list(all = mm, candidates = candidates[order(candidates$llik, decreasing = TRUE), 
        ])
}
.calcLlikSegmentExonLrs <- function(exon.lrs, log.ratio.offset, max.exon.ratio, sd.seg, 
    dt, b, D, test.num.copy) {
    c(.calcLlikSegmentSubClonal(exon.lrs + log.ratio.offset, max.exon.ratio), sapply(test.num.copy, 
        function(Ci) sum(dnorm(exon.lrs + log.ratio.offset, mean = log2(dt * Ci + 
            b/D), sd = sd.seg, log = TRUE))))
}

.rankResults <- function(results) {
    
    # max.ll <- max(sapply(results, function(z) z$log.likelihood)) max.snv.ll <-
    # max(sapply(results, function(z) z$SNV.posterior$beta.model$llik))
    
    for (i in 1:length(results)) {
        if (is.null(results[[i]]$SNV.posterior)) {
            results[[i]]$total.log.likelihood <- results[[i]]$log.likelihood
        } else {
            # results[[i]]$total.log.likelihood <- (results[[i]]$log.likelihood/max.ll) +
            # 2*(results[[i]]$SNV.posterior$beta.model$llik/max.snv.ll)
            results[[i]]$total.log.likelihood <- results[[i]]$log.likelihood + results[[i]]$SNV.posterior$beta.model$llik
        }
    }
    idx.opt <- order(sapply(results, function(z) z$total.log.likelihood), decreasing = TRUE)
    
    results <- results[idx.opt]
    
    # remove solutions with -inf likelihood score
    results <- results[!is.infinite(sapply(results, function(z) z$total.log.likelihood))]
    
    results
}
.sampleOffsetFast <- function(test.num.copy, seg, exon.lrs, sd.seg, p, C, total.ploidy, 
    max.exon.ratio, simulated.annealing, log.ratio.calibration = 0.25) {
    # Gibbs sample offset
    test.offset <- seq(sd.seg * -log.ratio.calibration, sd.seg * log.ratio.calibration, 
        by = 0.01)
    # adjust by chromosome or total. total works better in simulation? seg.ids.by.chr
    # <- lapply(split(cbind(seg, id=1:nrow(seg)), seg$chrom), function(x) x$id)
    seg.ids.by.chr <- list(1:nrow(seg))
    
    lr <- lapply(1:length(seg.ids.by.chr), function(j) {
        px.offset <- lapply(test.offset, function(px) sapply(seg.ids.by.chr[[j]], 
            function(i) {
                b <- 2 * (1 - p)
                D <- total.ploidy
                dt <- p/D
                llik.all <- .calcLlikSegmentExonLrs(exon.lrs[[i]], px, max.exon.ratio, 
                  sd.seg, dt, b, D, test.num.copy)
                sapply(llik.all, max)
            }))
        px.offset.s <- sapply(lapply(px.offset, apply, 2, max), sum)
        
        px.offset.s <- exp(px.offset.s - max(px.offset.s))
        log.ratio.offset <- test.offset[min(which(runif(n = 1, min = 0, max = sum(px.offset.s)) <= 
            cumsum(px.offset.s)))]
    })
    do.call(c, lapply(1:length(lr), function(i) rep(lr[[i]], length(seg.ids.by.chr[[i]]))))
}
.sampleOffset <- function(subclonal, seg, exon.lrs, sd.seg, p, C, total.ploidy, max.exon.ratio, 
    simulated.annealing, iter, log.ratio.calibration = 0.25) {
    # Gibbs sample offset
    test.offset <- seq(sd.seg * -log.ratio.calibration, sd.seg * log.ratio.calibration, 
        by = 0.01)
    # adjust by chromosome or total. total works better in simulation? seg.ids.by.chr
    # <- lapply(split(cbind(seg, id=1:nrow(seg)), seg$chrom), function(x) x$id)
    seg.ids.by.chr <- list(1:nrow(seg))
    
    lr <- lapply(1:length(seg.ids.by.chr), function(j) {
        px.offset <- lapply(test.offset, function(px) sapply(seg.ids.by.chr[[j]], 
            function(i) .calcLlikSegment(subclonal[i], exon.lrs[[i]] + px, sd.seg, 
                p, C[i], total.ploidy, max.exon.ratio)))
        
        px.offset.s <- sapply(px.offset, sum, na.rm = TRUE)
        if (simulated.annealing) 
            px.offset.s <- px.offset.s * exp(iter/4)
        px.offset.s <- exp(px.offset.s - max(px.offset.s))
        log.ratio.offset <- test.offset[min(which(runif(n = 1, min = 0, max = sum(px.offset.s)) <= 
            cumsum(px.offset.s)))]
    })
    do.call(c, lapply(1:length(lr), function(i) rep(lr[[i]], length(seg.ids.by.chr[[i]]))))
}
.smoothOutliers <- function(x, na.rm = TRUE, ...) {
    if (length(x) < 5) 
        return(x)
    qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm, ...)
    H <- 1.5 * IQR(x, na.rm = na.rm)
    y <- x
    # find points outside the 'boxplot' range
    idx <- x <= (qnt[1] - H) | x >= (qnt[2] + H)
    idx[c(1, length(idx))] <- FALSE
    # smooth outliers by using the average log-ratio including the two neighbors
    y[idx] <- (x[which(idx) - 1] + x[which(idx)] + x[which(idx) + 1])/3
    # are the first or last elements outliers? Then remove those (only one neighbor)
    idx <- x <= (qnt[1] - H) | x >= (qnt[2] + H)
    if (idx[1]) 
        y[1] <- NA
    if (idx[length(y)]) 
        y[length(y)] <- NA
    y[!is.na(y)]
}
.getGeneCallsLOH <- function(result) {
    g <- result$gene.calls
    g.gr <- GRanges(seqnames = g$chr, IRanges(start = g$start, end = g$end))
    l <- result$SNV.posterior$beta.model$loh$output
    l.gr <- GRanges(seqnames = l$chrom, IRanges(start = l$loc.start, end = l$loc.end))
    ov <- findOverlaps(g.gr, l.gr)
    g$num.snps.loh.segment <- 0
    g$percentage.loh.in.loh.segment <- 0
    g$num.snps.loh.segment[queryHits(ov)] <- l[subjectHits(ov), "num.mark"]
    g$percentage.loh.in.loh.segment[queryHits(ov)] <- l[subjectHits(ov), "seg.mean"]
    g
}
.getVariantCallsLOH <- function(result) {
    g <- result$SNV.posterior$beta.model$posteriors
    g.gr <- GRanges(seqnames = g$seqnames, IRanges(start = g$start, end = g$end))
    l <- result$SNV.posterior$beta.model$loh$output
    l.gr <- GRanges(seqnames = l$chrom, IRanges(start = l$loc.start, end = l$loc.end))
    ov <- findOverlaps(g.gr, l.gr)
    g$num.snps.loh.segment <- 0
    g$percentage.loh.in.loh.segment <- 0
    g$num.snps.loh.segment[queryHits(ov)] <- l[subjectHits(ov), "num.mark"]
    g$percentage.loh.in.loh.segment[queryHits(ov)] <- l[subjectHits(ov), "seg.mean"]
    
    g
}

.calcPuritySomaticVariants <- function(vcf, prior.somatic, tumor.id.in.vcf) {
    median(unlist(geno(vcf[prior.somatic > 0.5])$FA[, tumor.id.in.vcf]), na.rm = TRUE)/0.48
}
.createFakeLogRatios <- function(tumor, seg.file) {
    seg <- read.delim(seg.file)
    required.colnames <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", 
        "seg.mean")
    if (!all.equal(colnames(seg), required.colnames)) {
        stop(paste("Segmentation file expected with colnames", 
                paste(required.colnames, collapse = ", ")))
    }
    
    seg.gr <- GRanges(seqnames = .add.chr.name(seg$chrom), 
                IRanges(start = round(seg$loc.start), end = seg$loc.end))

    exon.gr <- GRanges(seqnames = gsub("24", "Y", gsub("23", "X", tumor$chr)), 
                IRanges(start = tumor$probe_start, end = tumor$probe_end))

    ov <- findOverlaps(exon.gr, seg.gr)
    log.ratio <- seg$seg.mean[subjectHits(ov)]
    # sanity check, so that every exon has exactly one segment log-ratio
    log.ratio <- log.ratio[match(1:nrow(tumor), queryHits(ov))]
    log.ratio
}
.strip.chr.name <- function(ls) {
    chr.hash <- NULL
    data(chr.hash, envir = environment())
    x <- chr.hash[as.character(ls), 2]
    x[is.na(x)] <- as.numeric(ls[is.na(x)])
    x
}
.add.chr.name <- function(ls) {
    chr.hash <- NULL
    data(chr.hash, envir = environment())
    x <- as.character(chr.hash$chr[match(ls, chr.hash$number)])
    x[is.na(x)] <- ls[is.na(x)]
    x
}
    
 
