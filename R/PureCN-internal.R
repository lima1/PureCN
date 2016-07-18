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
    sum(dunif(
        x =  vapply(2^lr, function(y) min(y, max.exon.ratio), double(1)), 
        min = 0, max = max.exon.ratio, log = TRUE))
}
.calcLogRatio <- function(normal, tumor, verbose) {
    # make sure that normal and tumor align
    if (!identical(as.character(normal[, 1]), as.character(tumor[, 1]))) {
        .stopUserError("Interval files in normal and tumor different.")
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

.calcMsSegmentC <- function(yy, test.num.copy, Ci, prior.K) {
    max.M <- floor(Ci/2)
    idx.germline <- test.num.copy+length(test.num.copy)+1
    idx.somatic <- test.num.copy+1

    yys <- lapply(0:max.M, function(Mi) {
        for (i in test.num.copy) {
            n.cases.germ <- ifelse(Mi==Ci-Mi,1,2)
            n.cases.somatic <- length(unique(c(1,Mi, Ci-Mi)))
             
            if (i!=Mi && i!=Ci - Mi) {
                yy[,idx.germline[i+1]] <- yy[,idx.germline[i+1]] + log(1-prior.K) - log(n.cases.germ)

                # allow somatic mutations always have M=1
                if (i==1) {
                    yy[,idx.somatic[i+1]] <- yy[,idx.somatic[i+1]] + log(prior.K) - log(n.cases.somatic)
                } else {
                    yy[,idx.somatic[i+1]] <- yy[,idx.somatic[i+1]] +
                        log(1-prior.K) - log(n.cases.somatic)
                }    
            } else {
                yy[,idx.germline[i+1]] <- yy[,idx.germline[i+1]] + log(prior.K) -log(n.cases.germ)
                yy[,idx.somatic[i+1]] <- yy[,idx.somatic[i+1]] + log(prior.K) - log(n.cases.somatic)
            }    
        }
        yy
    })    
    best <- which.max(vapply(yys, function(x) sum(apply(x, 1, max)),double(1)))
    list(best=yys[[best]], M=test.num.copy[best])
}

.calcMsSegment <- function(xxi, test.num.copy, prior.K) {
    lapply(seq_along(xxi), function(i).calcMsSegmentC( xxi[[i]], test.num.copy,
test.num.copy[i], prior.K))
}    

.calcSNVLLik <- function(vcf, tumor.id.in.vcf, ov, p, test.num.copy, 
    C.posterior, C, snv.model, prior.somatic, snv.lr, sampleid = NULL, 
    cont.rate = 0.01, prior.K) {

    prior.cont <- ifelse(info(vcf)$DB, cont.rate, 0)
    prior.somatic <- prior.somatic/(1 + cont.rate)
    
    haploid.penalty <- 0
    
    if (median(C) < 1.1 && p <= 0.3) {
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

    # the same issue, but for subclonal alterations below < 7
    idx <- C.posterior[,ncol(C.posterior)] > 0.999 & C < max(test.num.copy)
    if (length(idx) > 1) {
        f <- 1/ncol(C.posterior)
        C.posterior[idx,] <-(C.posterior[idx,]+f)/(ncol(C.posterior)*f+1)
    }   
    
    seg.idx <- which(seq_len(nrow(C.posterior)) %in% queryHits(ov))
    sd.ar <- sd(unlist(geno(vcf)$FA[, tumor.id.in.vcf]))
    xx <- lapply(seg.idx, function(i) {
        # classify germline vs somatic
        idx <- subjectHits(ov)[queryHits(ov) == i]
        
        ar_all <- unlist(geno(vcf)$FA[idx, tumor.id.in.vcf])
        ar_all[ar_all > 1] <- 1
        dp_all <- unlist(geno(vcf)$DP[idx, tumor.id.in.vcf])
        mInf_all <- log(double(length(ar_all)))
        shape1 <- ar_all * dp_all + 1
        shape2 <- (1 - ar_all) * dp_all + 1

        lapply(test.num.copy, function(Ci) {
            priorM <- log(Ci + 1 + haploid.penalty)
            
            skip <- test.num.copy > Ci | C.posterior[i, Ci + 1] <=0

            p.ar <- lapply(c(0, 1), function(g) {
                dbetax <- (p * test.num.copy + g * (1-p))/(p * Ci + 2 * (1-p))
                l <- lapply(seq_along(test.num.copy), function(j) {
                    if (skip[j]) return(mInf_all)
                    dbeta(x = dbetax[j], 
                    shape1 = shape1,  
                    shape2 = shape2, log = TRUE) - priorM
                })
                do.call(cbind,l)
            })
            
            p.ar.cont.1 <- dbeta(x = (p * Ci + 2 * (1 - p - cont.rate))/
                  (p * Ci + 2 * (1 - p)),shape1=shape1, shape2=shape2,log=TRUE) - 
                                  priorM

            p.ar.cont.2 <- dbeta(x = (cont.rate)/(p * Ci + 2 * (1 - p)), 
                  shape1=shape1, shape2=shape2,log=TRUE) -
                                  priorM

            # add prior probabilities for somatic vs germline
            p.ar[[1]] <- p.ar[[1]] + log(prior.somatic[idx])

            p.ar[[2]] <- p.ar[[2]] + log(1 - prior.somatic[idx])

            # contamination (either homozygous germline, or germline from 
            # other sample)

            p.ar[[3]] <- p.ar.cont.1 + log(prior.cont[idx])
            p.ar[[4]] <- p.ar.cont.2 + log(prior.cont[idx])
            
            do.call(cbind, p.ar)
        })
        
    })
    
    tmp <- lapply(xx, .calcMsSegment, test.num.copy, prior.K)
    xx <- lapply(tmp, lapply, function(x) x$best)
    
    # Get segment M's for each SNV
    segment.M <- sapply(seq_along(tmp), function(i) 
        tmp[[i]][[min(C[seg.idx[i]], max(test.num.copy))+1]]$M)
    segment.M <- unlist(sapply(seq_along(seg.idx), function(i) 
        rep(segment.M[i], sum(seg.idx[i]==queryHits(ov)))))

    snv.posteriors <- do.call(rbind, 
        lapply(seq_along(xx), function(i) Reduce("+", 
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
        ML.C = C[queryHits(ov)],
        ML.M.Segment=segment.M
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
            presorted = TRUE),verbose=0)
    
    # these are potential artifacts with very high clonal probability and would have
    # huge impact on log-likelihood
    rm.snv.posteriors <- apply(snv.posteriors, 1, max)
    idx.ignore <- (posteriors$CN.Subclonal & 
        posteriors$ML.C < max(test.num.copy) & 
        C.posterior[queryHits(ov), ncol(C.posterior)] > 0.95) | 
        rm.snv.posteriors == 0

    # change seqnames to chr
    colnames(posteriors)[1] <- "chr"    

    ret <- list(
        llik = sum(log(rm.snv.posteriors[!idx.ignore])) - sum(idx.ignore), 
        likelihoods = snv.posteriors, 
        posteriors = posteriors, 
        vcf.ids = vcf.ids, 
        segment.ids = queryHits(ov), 
        loh = loh, 
        llik.ignored = idx.ignore)

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

.checkParameters <- function(test.purity, min.ploidy, max.ploidy, 
    max.non.clonal) {
    if (min(test.purity) <= 0 || max(test.purity) > 1) 
        .stopUserError("test.purity not within expected range.")
    if (min.ploidy <= 0 || max.ploidy <= 2) 
        .stopUserError("min.ploidy or max.ploidy not within expected range.")
    if (max.non.clonal > 1) 
        .stopUserError("max.non.clonal not within expected range.")

    stopifnot(is.numeric(min.ploidy))
    stopifnot(is.numeric(max.ploidy))
    stopifnot(is.numeric(test.purity))
    stopifnot(is.numeric(max.non.clonal))
}
.failedNonAberrant <- function(result, cutoffs = c(0.01, 0.005)) {
    xx <- split(result$seg, result$seg$C)
    if (length(xx) < 3) 
        return(TRUE)
    xx.sum <- sort(vapply(xx, function(x) sum(x$size), double(1)), 
        decreasing = TRUE)
    xx.sum <- xx.sum/sum(xx.sum)
    if (xx.sum[2] <= cutoffs[1] && xx.sum[3] <= cutoffs[2]) 
        return(TRUE)
    FALSE
}
.filterDuplicatedResults <- function(results) {
    if (length(results) < 2) 
        return(results)
    idx.duplicated <- rep(FALSE, length(results))

    for (i in seq_len(length(results)-1)) {
        for (j in seq(i+1,length(results)) ) {
            if ( abs( results[[i]]$purity - results[[j]]$purity ) < 0.1 &&
                 abs( results[[i]]$ploidy - results[[j]]$ploidy ) / results[[i]]$ploidy < 0.1) {
                idx.duplicated[j] <- TRUE
            }    
        }    
    }
    results[!idx.duplicated]
}
.findLocalMinima <- function(m) {
    loc.min <- matrix(nrow = 0, ncol = 2)
    for (i in seq_len(nrow(m))) {
        for (j in seq_len(ncol(m))) {
            x <- seq(i - 1, i + 1)
            x <- x[x >= 1 & x <= nrow(m)]
            y <- seq(j - 1, j + 1)
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
    fraction.loh <- .getFractionLoh(result)
    # Assume that everything below 2.6 did not undergo genome duplication, which can
    # result in lots of LOH
    if (result$ploidy < 2.6 && fraction.loh > 0.5) {
        result$flag <- TRUE
        result$flag_comment <- .appendComment(result$flag_comment, "EXCESSIVE LOH")
    }
    return(result)
}

.getFractionLoh <- function(result, min.seg.mean=0.3) {
    tmp <- result$SNV.posterior$beta.model$loh$output
    if (is.null(tmp)) return(0)
    tmp <- tmp[complete.cases(tmp), ]
    sum(tmp[which(tmp$seg.mean >= min.seg.mean), "num.mark"])/sum(tmp$num.mark)
}

.flagResults <- function(results, max.non.clonal = 0.2, max.logr.sdev, logr.sdev, max.segments,
    flag = NA, flag_comment = NA, dropout=FALSE) {
    results <- lapply(results, .flagResult, max.non.clonal = max.non.clonal)
    
    # ldiff <- ( results[[1]]$total.log.likelihood -
    # results[[2]]$total.log.likelihood ) / abs( results[[1]]$total.log.likelihood )
    # if (ldiff < 0.1) { results[[1]]$flag <- TRUE results[[1]]$flag_comment <-
    # 'SIMILAR TO SECOND MOST LIKELY OPTIMUM' return(results) }
    number.segments <- nrow(results[[1]]$seg)
    
    if (logr.sdev > max.logr.sdev) {
        for (i in seq_along(results)) {
            results[[i]]$flag <- TRUE
            results[[i]]$flag_comment <- .appendComment(results[[i]]$flag_comment, 
                "NOISY LOG-RATIO")
        }
    }

    if (number.segments > max.segments) {
        for (i in seq_along(results)) {
            results[[i]]$flag <- TRUE
            results[[i]]$flag_comment <- .appendComment(results[[i]]$flag_comment, 
                "NOISY SEGMENTATION")
        }
    }

    if (dropout) {
        for (i in seq_along(results)) {
            results[[i]]$flag <- TRUE
            results[[i]]$flag_comment <- .appendComment(results[[i]]$flag_comment, 
                "HIGH AT- OR GC-DROPOUT")
        }
    }
    
    # some global flags
    if (!is.na(flag) && flag) {
        for (i in seq_along(results)) {
            results[[i]]$flag <- TRUE
            results[[i]]$flag_comment <- .appendComment(results[[i]]$flag_comment, 
                flag_comment)
        }
    }
    results
}
.getGeneCalls <- function(seg.adjusted, gc.data, log.ratio, fun.focal, args.focal, chr.hash) {
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
    
    gc.gr <- GRanges(seqnames = .strip.chr.name(gc.pos$chrom, chr.hash), IRanges(start = gc.pos$start, 
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
            llik.all <- lapply(seq_along(exon.lrs), function(i) .calcLlikSegmentExonLrs(exon.lrs[[i]], 
                log.ratio.offset, max.exon.ratio, sd.seg, dt, b, D, test.num.copy))
            subclonal <- vapply(llik.all, which.max, double(1)) == 1
            subclonal.f <- length(unlist(exon.lrs[subclonal]))/length(unlist(exon.lrs))
            if (debug) 
                message(paste(sum(subclonal), subclonal.f))
            if (subclonal.f > max.non.clonal) 
                return(-Inf)
            llik <- sum(vapply(llik.all, max, double(1)))
        })
    })
    mm <- sapply(mm, function(x) unlist(x))
    colnames(mm) <- test.purity
    rownames(mm) <- ploidy.grid
    ai <- .findLocalMinima(mm)
    candidates <- data.frame(ploidy = as.numeric(rownames(mm)[ai[, 1]]), purity = as.numeric(colnames(mm)[ai[, 
        2]]), llik = mm[ai])
    candidates$tumor.ploidy <- (candidates$ploidy - 2 * (1 - candidates$purity))/candidates$purity
    
    # add diploid candidate solutions in the following purity grid in 
    # case there are none.
    grid <- seq(0,1,by=1/4)
    for (i in seq_along(grid)[-length(grid)]) {
   
        t1 <- which.min(abs(as.numeric(colnames(mm)) - grid[i]))
        t2 <- which.min(abs(as.numeric(colnames(mm)) - grid[i+1]))
        if (t2-t1 < 2) next

        # Nothing close to diplod in this range? Then add.
        if (min(abs(2 - candidates$tumor.ploidy[candidates$purity>grid[i] & 
            candidates$purity< grid[i+1] ])) > 0.3) {

            mm.05 <- mm[,(t1+1):t2,drop=FALSE]
            candidates <- rbind(candidates, 
                c(2, as.numeric(names(which.max(mm.05["2", ]))), 
                max(mm.05["2", ]), 2))

            # Remove again if too similar with existing candidate
            if (nrow(candidates) > 2 && 
                abs(Reduce("-",tail(candidates$ploidy,2))) < 0.001 && 
                abs(Reduce("-",tail(candidates$purity,2))) < 0.1) {
                candidates <- candidates[- (nrow(candidates) - 2 + 
                    which.min(tail(candidates$llik,2))),]
            }    
        }
    }
    
    candidates <- candidates[candidates$tumor.ploidy >= 0.5, ]
    
    list(all = mm, candidates = candidates[order(candidates$llik, decreasing = TRUE), 
        ])
}
.calcLlikSegmentExonLrs <- function(exon.lrs, log.ratio.offset, max.exon.ratio, sd.seg, 
    dt, b, D, test.num.copy) {
    c(.calcLlikSegmentSubClonal(exon.lrs + log.ratio.offset, max.exon.ratio), vapply(test.num.copy, 
        function(Ci) sum(dnorm(exon.lrs + log.ratio.offset, mean = log2(dt * Ci + 
            b/D), sd = sd.seg, log = TRUE)),double(1)))
}

# This function is used to punish more complex copy number models
# a little bit. Based on the BIC. This just counts the number of utilized 
# copy number states, excluding normal 2. Then multiplies by 
# log(number exons)
.calcComplexityCopyNumber <- function(results) {
    cs <- sapply((0:7)[-3], function(i) sapply(results, function(y)
                    sum(y$seg$size[y$seg$C == i])/sum(y$seg$size)))
    complexity <- apply(cs,1, function(z) sum(z>0.001))
    n <- sum(results[[1]]$seg$num.mark, na.rm=TRUE)
    penalty <- -complexity*log(n)
}

.rankResults <- function(results) {
    
    # max.ll <- max(sapply(results, function(z) z$log.likelihood)) max.snv.ll <-
    # max(sapply(results, function(z) z$SNV.posterior$beta.model$llik))
    complexity <- .calcComplexityCopyNumber(results) 
    for (i in seq_along(results)) {
        if (is.null(results[[i]]$SNV.posterior)) {
            results[[i]]$total.log.likelihood <- results[[i]]$log.likelihood
        } else {
            # results[[i]]$total.log.likelihood <- (results[[i]]$log.likelihood/max.ll) +
            # 2*(results[[i]]$SNV.posterior$beta.model$llik/max.snv.ll)
            results[[i]]$total.log.likelihood <- results[[i]]$log.likelihood + results[[i]]$SNV.posterior$beta.model$llik
        }
        results[[i]]$total.log.likelihood <- results[[i]]$total.log.likelihood + complexity[i]
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
    seg.ids.by.chr <- list(seq_len(nrow(seg)))
    
    lr <- lapply(seq_along(seg.ids.by.chr), function(j) {
        px.offset <- lapply(test.offset, function(px) vapply(seg.ids.by.chr[[j]], 
            function(i) {
                b <- 2 * (1 - p)
                D <- total.ploidy
                dt <- p/D
                llik.all <- .calcLlikSegmentExonLrs(exon.lrs[[i]], px, max.exon.ratio, 
                  sd.seg, dt, b, D, test.num.copy)
                vapply(llik.all, max, double(1))
            }, double(1+length(test.num.copy))))
        px.offset.s <- vapply(lapply(px.offset, apply, 2, max), sum, double(1))
        
        px.offset.s <- exp(px.offset.s - max(px.offset.s))
        log.ratio.offset <- test.offset[min(which(runif(n = 1, min = 0, max = sum(px.offset.s)) <= 
            cumsum(px.offset.s)))]
    })
    do.call(c, lapply(seq_along(lr), function(i) rep(lr[[i]], length(seg.ids.by.chr[[i]]))))
}
.sampleOffset <- function(subclonal, seg, exon.lrs, sd.seg, p, C, total.ploidy, max.exon.ratio, 
    simulated.annealing, iter, log.ratio.calibration = 0.25) {
    # Gibbs sample offset
    test.offset <- seq(sd.seg * -log.ratio.calibration, sd.seg * log.ratio.calibration, 
        by = 0.01)
    seg.ids.by.chr <- list(seq_len(nrow(seg)))
    
    lr <- lapply(seq_along(seg.ids.by.chr), function(j) {
        px.offset <- lapply(test.offset, function(px) vapply(seg.ids.by.chr[[j]], 
            function(i) .calcLlikSegment(subclonal[i], exon.lrs[[i]] + px, sd.seg, 
                p, C[i], total.ploidy, max.exon.ratio), double(1))
        )
        
        px.offset.s <- sapply(px.offset, sum, na.rm = TRUE)
        if (simulated.annealing) 
            px.offset.s <- px.offset.s * exp(iter/4)
        px.offset.s <- exp(px.offset.s - max(px.offset.s))
        log.ratio.offset <- test.offset[min(which(runif(n = 1, min = 0, max = sum(px.offset.s)) <= 
            cumsum(px.offset.s)))]
    })
    do.call(c, lapply(seq_along(lr), function(i) rep(lr[[i]], length(seg.ids.by.chr[[i]]))))
}
.removeOutliers <- function(x, na.rm=TRUE,...) {
    if (length(x) < 5) 
        return(x)
    qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm, ...)
    H <- 1.5 * IQR(x, na.rm = na.rm)
    if (is.na(H) || H < 0.001) return(x)
    # find points outside the 'boxplot' range
    idx <- x <= (qnt[1] - H) | x >= (qnt[2] + H)
    x[!idx]
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

.calcPuritySomaticVariants <- function(vcf, prior.somatic, tumor.id.in.vcf) {
    median(unlist(geno(vcf[prior.somatic > 0.5])$FA[, tumor.id.in.vcf]), na.rm = TRUE)/0.48
}
.createFakeLogRatios <- function(tumor, seg.file, chr.hash) {
    seg <- read.delim(seg.file)
    required.colnames <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", 
        "seg.mean")
    if (!identical(colnames(seg), required.colnames)) {
        .stopUserError(paste("Segmentation file expected with colnames", 
                paste(required.colnames, collapse = ", ")))
    }
    
    seg.gr <- GRanges(seqnames = .add.chr.name(seg$chrom, chr.hash), 
                IRanges(start = round(seg$loc.start), end = seg$loc.end))

    exon.gr <- GRanges(seqnames = tumor$chr, 
                IRanges(start = tumor$probe_start, end = tumor$probe_end))

    ov <- findOverlaps(exon.gr, seg.gr)
    log.ratio <- seg$seg.mean[subjectHits(ov)]
    # sanity check, so that every exon has exactly one segment log-ratio
    log.ratio <- log.ratio[match(seq_len(nrow(tumor)), queryHits(ov))]
    mean.log.ratio <- mean(subset(log.ratio, !is.infinite(log.ratio)), 
        na.rm = TRUE)
    # calibrate
    log.ratio <- log.ratio - mean.log.ratio
    log.ratio
}
.strip.chr.name <- function(ls, chr.hash) {
    x <- chr.hash[as.character(ls), 2]
    x[is.na(x)] <- as.numeric(ls[is.na(x)])
    x
}
.add.chr.name <- function(ls, chr.hash) {
    x <- as.character(chr.hash$chr[match(ls, chr.hash$number)])
    x[is.na(x)] <- ls[is.na(x)]
    x
}
.getChrHash <- function(ls) {
    ls <- unique(ls)
    ls <- ls[!ls %in% "chrM"]
    chr.hash <- NULL
    data(chr.hash, envir = environment())

    if (sum(!ls %in% chr.hash[,1]) == 0) return(chr.hash)

    data.frame(chr=as.factor(ls), number=seq_along(ls), row.names=ls)
}
.checkChrHash <- function(chr.hash) {
    if (!identical(colnames(chr.hash)[1:2], c("chr", "number"))) {
        .stopUserError("Colnames of chr.hash should be 'chr' and 'number'.")
    }
    if (sum(duplicated(chr.hash$chr))>0) {
        .stopUserError("Duplicate chromosome names in chr.hash.")
    }
    rownames(chr.hash) <- chr.hash$chr
    chr.hash
}
#.ffpeCleanLogRatio <- function(log.ratio, window=20) {
#   dlr <- c(0, diff(log.ratio))
#   start <- seq(1,length(dlr), by=window)
#   end <- seq(window,length(dlr), by=window)
#   sds <- sapply(seq_along(end), function(i) sd(dlr[start[i]:end[i]], na.rm=TRUE))
##   ids <- which(sds > quantile(sds, na.rm=TRUE, p=1-0.001))
#}

# calculate allelic fraction from read depths
.addFaField <- function(vcf, field="FA") {
    if (is.null(geno(vcf)$AD)) {
        .stopRuntimeError("No AD geno field in VCF.")
    }
    matrixFA <- do.call(rbind, apply(geno(vcf)$AD,1, function(x) lapply(x,
        function(y) y[2]/sum(y))))
    newGeno <- DataFrame(
        Number="A", Type="Float",
        Description="Allele fraction of the alternate allele",
        row.names=field)
    geno(header(vcf)) <- rbind(geno(header(vcf)), newGeno)
    geno(vcf)[[field]] <- matrixFA
    vcf
}
.addDpField <- function(vcf, field="DP") {
    if (is.null(geno(vcf)$AD)) {
        .stopRuntimeError("No AD geno field in VCF.")
    }
    matrixDP <- apply(geno(vcf)$AD,2,sapply, sum)
    newGeno <- DataFrame(
        Number=1, Type="Integer",
        Description="Approximate read depth",
        row.names=field)
    geno(header(vcf)) <- rbind(geno(header(vcf)), newGeno)
    geno(vcf)[[field]] <- matrixDP
    vcf
}
.readAndCheckVcf <- function(vcf.file, genome) {
    if (class(vcf.file) == "character") {    
        vcf <- readVcf(vcf.file, genome)
    } else if (class(vcf.file) != "CollapsedVCF") {
        .stopUserError("vcf.file neither a filename nor a CollapsedVCF ", 
            "object.") 
    } else {
        vcf <- vcf.file
    } 
    if (is.null(info(vcf)$DB)) {
        .stopUserError(vcf.file, " has no DB info field for dbSNP membership.")
    }
    if (is.null(geno(vcf)$AD)) {
        .stopUserError(vcf.file, 
            " has no AD geno field containing read depths of ref and alt.")
    }
    if (is.null(geno(vcf)$FA)) {
        # try to add an FA geno field if missing
        vcf <- .addFaField(vcf)
    }
    if (is.null(geno(vcf)$DP)) {
        # try to add an FA geno field if missing
        vcf <- .addDpField(vcf)
    }
    vcf     
}    

.stopUserError <- function(...) {
    msg <- paste(c(...), collapse="")
    msg <- paste(msg, "\n\nThis is most likely a user error due to",
        " invalid input data or parameters (PureCN ", 
        packageVersion("PureCN"), ").", sep="")
    stop( paste(strwrap(msg),"\n"), call.= FALSE)
}
.stopRuntimeError <- function(...) {
    msg <- paste(c(...), collapse="")
    msg <- paste(msg, "\n\nThis runtime error might be caused by",
        " invalid input data or parameters. Please report bug (PureCN ", 
        packageVersion("PureCN"), ").", sep="")
    stop( paste(strwrap(msg),"\n"), call.= FALSE)
}

.calcGCmetric <- function(gc.data, coverage) { 
    gcbins <- split(coverage$average.coverage, gc.data$gc_bias < 0.5); 
    mean(gcbins[[1]], na.rm=TRUE)/mean(gcbins[[2]], na.rm=TRUE) 
}
.checkGCBias <- function(normal, tumor, gc.data, max.dropout, verbose) {
    gcMetricNormal <- .calcGCmetric(gc.data, normal)
    gcMetricTumor <- .calcGCmetric(gc.data, tumor)
    if (verbose) { 
        message("AT/GC dropout: ", round(gcMetricTumor, digits=2),
        " (tumor), ", round(gcMetricNormal, digits=2), " (normal).")
    }    
    if (gcMetricNormal < max.dropout[1] || 
        gcMetricNormal > max.dropout[2] ||
        gcMetricTumor  < max.dropout[1] ||
        gcMetricTumor  > max.dropout[2]) {
        warning("High GC-bias in normal or tumor. Is data GC-normalized?")
        return(TRUE)
    }
    return(FALSE)
}

.gcGeneToCoverage <- function(gc.gene.file, coverage.cutoff) {
    gc.data <- readCoverageGatk(gc.gene.file)
    gc.data$average.coverage <- coverage.cutoff
    gc.data$coverage <- coverage.cutoff * gc.data$targeted.base
    gc.data
}
.getSex <- function(sex, normal, tumor) {
    if (sex != "?") return(sex)
    sex.tumor <- getSexFromCoverage(tumor, verbose=FALSE)
    sex.normal <- getSexFromCoverage(normal, verbose=FALSE)
    if (!identical(sex.tumor, sex.normal)) {
        warning("Sex tumor/normal mismatch: tumor = ", sex.tumor, 
            " normal = ", sex.normal)
    }
    sex <- sex.tumor    
    if (is.na(sex)) sex = "?"
    sex
}
.fixAllosomeCoverage <- function(sex, tumor) {
    sex.chr <- .getSexChr(tumor$chr)
    if (sex=="M" || sex=="?") {
        tumor <- .removeChr(tumor, remove.chrs=sex.chr)
    } else if (sex=="F") {
        tumor <- .removeChr(tumor, remove.chrs=sex.chr[2])
    }       
    tumor
}
.filterExonsChrHash <- function(exonsUsed, tumor, chr.hash, verbose) {
    if (is.null(chr.hash)) return(exonsUsed)
    nBefore <- sum(exonsUsed)
    exonsUsed <- exonsUsed & tumor$chr %in% chr.hash$chr
    nAfter <- sum(exonsUsed)

    if (verbose && nAfter < nBefore) { 
        message("Removing ", nBefore-nAfter, " exons on chromosomes ",
            "outside chr.hash.")
    }
    exonsUsed
}
.filterExonsTargetedBase <- function(exonsUsed, tumor, filter.targeted.base,
    verbose) {
    if (is.null(filter.targeted.base)) return(exonsUsed)
    nBefore <- sum(exonsUsed)
    exonsUsed <- exonsUsed & !is.na(tumor$targeted.base) & 
        tumor$targeted.base >= filter.targeted.base
    nAfter <- sum(exonsUsed)
    if (verbose && nAfter < nBefore) message("Removing ", nBefore-nAfter, 
        " small exons.")
    exonsUsed
}
.filterExonsLohHighGC <- function(exonsUsed, tumor, gc.data,
    filter.lowhigh.gc.exons, verbose) {
    gc.data <- gc.data[match(as.character(tumor[,1]), gc.data[,1]),]
    qq <- quantile(gc.data$gc_bias, p=c(filter.lowhigh.gc.exons, 
        1-filter.lowhigh.gc.exons), na.rm=TRUE)

    nBefore <- sum(exonsUsed)
    exonsUsed <- exonsUsed & !is.na(gc.data$gc_bias) & 
        !(gc.data$gc_bias < qq[1] | gc.data$gc_bias > qq[2])
    nAfter <- sum(exonsUsed)

    if (verbose && nAfter < nBefore) message("Removing ", 
        nBefore-nAfter, " low/high GC exons.")

    exonsUsed
}


