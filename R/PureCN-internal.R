# Make CMD check happy
globalVariables(names=c("Gene", "LR", "chrom", "seg.id", "seg.length",
"seg.mean"))

# calculates the log-likelihood of a segment, given log-ratios, the
# log ratio standard deviation, purity, copy number and ploidy.
# for sub-clonal alterations, it uses a uniform distribution, otherwise
# multiple gaussians for all tested copy numbers. 
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
# Previously calculated likelihood scores do not take the segmentation into
# account. This will find the most likely segment minor copy number
.calcMsSegmentC <- function(yy, test.num.copy, Ci, prior.K, mapping.bias.ok, 
    seg.id, min.variants.segment) {
    prior.M <- list(0.2,0.15,c(0.1,0.25),c(0.1,0.3),c(0.1,0.2,0.55))
    prior.M <- c(list(1), lapply(prior.M, function(x) c(x, 1-sum(x))))
    max.M <- floor(Ci/2)
    idx.germline <- test.num.copy+length(test.num.copy)+1
    idx.somatic <- test.num.copy+1
    yys <- lapply(0:max.M, function(Mi) {
        for (i in test.num.copy) {
            n.cases.germ <- ifelse(Mi==Ci-Mi,1,2)
            tmp <- unique(c(0, 1,Mi, Ci-Mi))
            n.cases.somatic <- length(tmp)
            Cx <- max(i,Ci)

            if (i!=Mi && i!=Ci - Mi) {
                yy[,idx.germline[i+1]] <- yy[,idx.germline[i+1]] + log(1-prior.K) - log(Cx + 1 - n.cases.germ)

                # allow somatic mutations always have M=1
                if (i<=1) {
                    yy[,idx.somatic[i+1]] <- yy[,idx.somatic[i+1]] + log(prior.K) - log(n.cases.somatic)
                } else {
                    yy[,idx.somatic[i+1]] <- yy[,idx.somatic[i+1]] +
                        log(1-prior.K) - log(Cx + 1 - n.cases.somatic)
                }    
            } else {
                yy[,idx.germline[i+1]] <- yy[,idx.germline[i+1]] + log(prior.K) -log(n.cases.germ)
                yy[,idx.somatic[i+1]] <- yy[,idx.somatic[i+1]] + log(prior.K) - log(n.cases.somatic)
            }    
        }
        yy
    })    
    # if not enough variants in segment, flag
    flag <- sum(mapping.bias.ok) < min.variants.segment
    if (!sum(mapping.bias.ok)) {
        mapping.bias.ok <- rep(TRUE, length(mapping.bias.ok))
    }

    likelihoodScores <- vapply(yys, function(x) 
        sum(apply(x[mapping.bias.ok,,drop=FALSE], 1, max)),double(1))
    best <- which.max(likelihoodScores)
    # transform and scale 
    likelihoodScores <- exp(likelihoodScores - likelihoodScores[best])
    posterior <- likelihoodScores[best]/sum(likelihoodScores, na.rm=TRUE)
    if (is.na(posterior) || posterior < 0.5) flag <- TRUE 

    list(best=yys[[best]], M=test.num.copy[best], flag=flag, posterior=posterior)
}

# calculate likelihood scores for all possible minor copy numbers
.calcMsSegment <- function(xxi, test.num.copy, opt.C, prior.K, mapping.bias.ok, 
    seg.id, min.variants.segment) {
    lapply(seq_along(xxi), function(i).calcMsSegmentC(xxi[[i]], test.num.copy,
c(test.num.copy, round(opt.C))[i], prior.K, mapping.bias.ok, seg.id, min.variants.segment))
}    

.calcSNVLLik <- function(vcf, tumor.id.in.vcf, ov, p, test.num.copy, 
    C.likelihood, C, opt.C, median.C, snv.model, prior.somatic, mapping.bias, snv.lr, 
    sampleid = NULL, cont.rate = 0.01, prior.K, max.coverage.vcf, non.clonal.M,
    model.homozygous=FALSE, error=0.001, max.mapping.bias=0.8, max.pon,
    min.variants.segment) {
    
    .dbeta <- function(x, shape1, shape2, log, size) dbeta(x=x, 
        shape1=shape1, shape2=shape2, log=log)
    if (snv.model=="betabin") {
        .dbeta <- function(x, shape1, shape2, log, size) {
            size <- pmin(size, max.coverage.vcf)
            dbetabinom.ab(x=round(x*size),shape1=shape1, shape2=shape2, 
                size=size, log=TRUE)
         }
    }

    prior.cont <- ifelse(prior.somatic < 0.1, cont.rate, 0)
    prior.somatic <- prior.somatic - (prior.cont*prior.somatic)
    priorHom <- if (model.homozygous) -log(3) else log(0)

    if (length(prior.somatic) != nrow(vcf)) {
        .stopRuntimeError("prior.somatic and VCF do not align.")
    }
    if (length(mapping.bias$bias) != nrow(vcf)) {
        .stopRuntimeError("mapping.bias and VCF do not align.")
    }
    
    haploid.penalty <- 0
    
    if (median.C < 1.1) {
        haploid.penalty <- 2
    }
    
    subclonal <- apply(C.likelihood[queryHits(ov), ], 1, which.max) == ncol(C.likelihood)
    
    seg.idx <- which(seq_len(nrow(C.likelihood)) %in% queryHits(ov))

    ar_all <- unlist(geno(vcf)$FA[, tumor.id.in.vcf])
    ar_all <- ar_all/mapping.bias$bias
    ar_all[ar_all > 1] <- 1
    
    dp_all <- unlist(geno(vcf)$DP[, tumor.id.in.vcf])
    if (snv.model!="betabin") dp_all[dp_all>max.coverage.vcf] <- max.coverage.vcf

    # Fit variants in all segments
    xx <- lapply(seg.idx, function(i) {
        # classify germline vs somatic
        idx <- subjectHits(ov)[queryHits(ov) == i]
        
        shape1 <- ar_all[idx] * dp_all[idx] + 1
        shape2 <- (1 - ar_all[idx]) * dp_all[idx] + 1
        mInf_all <- log(double(length(shape1)))

        list(vcf.ids=idx, likelihoods=lapply(seq(ncol(C.likelihood)), function(k) {
            Ci <- c(test.num.copy, opt.C[i])[k]
            priorM <- log(max(Ci,1) + 1 + haploid.penalty)
            
            skip <- test.num.copy > Ci | C.likelihood[i, k] <= 0

            p.ar <- lapply(c(0, 1), function(g) {
                cns <- test.num.copy
                if (!g) cns[1] <- non.clonal.M
                dbetax <- (p * cns + g * (1-p)) / (p * Ci + 2 * (1-p))
                l <- lapply(seq_along(test.num.copy), function(j) {
                    if (skip[j]) return(mInf_all)
                    .dbeta(x = dbetax[j],
                    shape1 = shape1,
                    shape2 = shape2, log = TRUE, size = dp_all[idx]) - priorM
                })
                do.call(cbind,l)
            })
            
            p.ar.cont.1 <- .dbeta(x = (p * Ci + 2 * (1 - p - cont.rate))/
                  (p * Ci + 2 * (1 - p)),shape1=shape1, shape2=shape2, 
                  log=TRUE, size=dp_all[idx]) - priorM

            p.ar.cont.2 <- .dbeta(x = cont.rate / (p * Ci + 2 * (1 - p)), 
                  shape1 = shape1, shape2 = shape2, log = TRUE, 
                  size = dp_all[idx]) - priorM

            # add prior probabilities for somatic vs germline
            p.ar[[1]] <- p.ar[[1]] + log(prior.somatic[idx])

            p.ar[[2]] <- p.ar[[2]] + log(1 - prior.somatic[idx])

            # contamination (either homozygous germline, or germline from 
            # other sample)

            p.ar[[3]] <- p.ar.cont.1 + log(prior.cont[idx])
            p.ar[[4]] <- p.ar.cont.2 + log(prior.cont[idx])

            # homozygous state
            p.ar[[5]] <- dbinom(round((1-ar_all[idx])*dp_all[idx]),size=round(dp_all[idx]), 
                prob=error/3, log=TRUE) + priorHom + log(1-prior.somatic[idx])
             
            do.call(cbind, p.ar)
        }))
        
    })
    
    tmp <- lapply(seq_along(xx),function(i) .calcMsSegment(xx[[i]]$likelihoods, 
               test.num.copy, opt.C[seg.idx[i]], prior.K, 
               mapping.bias.ok=mapping.bias$bias[xx[[i]]$vcf.ids]>=max.mapping.bias,
               seg.id=seg.idx[i], min.variants.segment))

    xx <- lapply(tmp, lapply, function(x) x$best)
    
    .extractValues <- function(tmp, field) {
        segmentValue <- sapply(seq_along(tmp), function(i) 
            tmp[[i]][[min(C[seg.idx[i]], max(test.num.copy))+1]][[field]])
        segmentValue <- unlist(sapply(seq_along(seg.idx), function(i) 
            rep(segmentValue[i], sum(seg.idx[i]==queryHits(ov)))))
    }
    # Get segment M's for each SNV
    segment.M <- .extractValues(tmp, "M")
    segment.Flag <- .extractValues(tmp, "flag")
    segment.Posterior <- .extractValues(tmp, "posterior")

    likelihoods <- do.call(rbind, 
        lapply(seq_along(xx), function(i) Reduce("+", 
            lapply(seq(ncol(C.likelihood)), function(j) 
                exp(xx[[i]][[j]]) * C.likelihood[seg.idx[i], j]))))

    colnames(likelihoods) <- c(paste("SOMATIC.M", test.num.copy, sep = ""), 
        paste("GERMLINE.M", test.num.copy, sep = ""), "GERMLINE.CONTHIGH", 
        "GERMLINE.CONTLOW", "GERMLINE.HOMOZYGOUS")
    
    vcf.ids <- do.call(c, lapply(seg.idx, function(i) 
        subjectHits(ov)[queryHits(ov) == i]))
    rownames(likelihoods) <- vcf.ids
    
    # for very high-level amplifications, all posteriors can be 0, so make sure
    # we get valid values here and flag those later.
    posteriors <- likelihoods/pmax(rowSums(likelihoods),.Machine$double.xmin)
    # this just adds a lot of helpful info to the SNV posteriors
    xx <- .extractMLSNVState(posteriors)
    
    posteriors <- cbind(
        as.data.frame(rowRanges(vcf[vcf.ids]), row.names=NULL)[, 1:3],
        ID = names(vcf[vcf.ids]),
        REF = as.character(ref(vcf[vcf.ids])),
        ALT = sapply(alt(vcf[vcf.ids]), function(x) 
                     paste(as.character(x), collapse=";")),
        posteriors,
        xx, 
        ML.C = C[queryHits(ov)],
        ML.M.SEGMENT = segment.M,
        M.SEGMENT.POSTERIOR = segment.Posterior,
        M.SEGMENT.FLAGGED = segment.Flag,
        row.names = NULL
    )
    
    posteriors$ML.AR <- (p * posteriors$ML.M + 
        ifelse(posteriors$ML.SOMATIC, 0, 1) * 
        (1 - p)) / (p * posteriors$ML.C + 2 * (1 - p))
    posteriors$ML.AR[posteriors$ML.AR > 1] <- 1 

    posteriors$AR <- unlist(geno(vcf[vcf.ids])$FA[, tumor.id.in.vcf])
    posteriors$AR.ADJUSTED <- posteriors$AR / mapping.bias$bias[vcf.ids]
    posteriors$AR.ADJUSTED[posteriors$AR.ADJUSTED>1] <- 1
    posteriors$MAPPING.BIAS <- mapping.bias$bias[vcf.ids]
    # Extract LOH
    posteriors$ML.LOH <- (posteriors$ML.M == posteriors$ML.C | 
        posteriors$ML.M == 0 | posteriors$ML.C == 1)
    
    posteriors$CN.SUBCLONAL <- subclonal
    depth <-as.numeric(geno(vcf[vcf.ids])$DP[, tumor.id.in.vcf])
    ar <- posteriors$AR.ADJUSTED
    ar[!posteriors$ML.SOMATIC] <- NA

    m <-  t(apply(cbind(ar, depth,  posteriors$ML.C), 1, function(x) 
        .calculate_ccf(vaf = x[1], depth = x[2], purity = p, C = x[3])))

    posteriors$CELLFRACTION <- as.numeric(m[,1])
    posteriors$CELLFRACTION.95.LOWER <- as.numeric(m[,2])
    posteriors$CELLFRACTION.95.UPPER <- as.numeric(m[,3])

    rm.snv.posteriors <- apply(likelihoods, 1, max)
    idx.ignore <- rm.snv.posteriors == 0 |
        posteriors$MAPPING.BIAS < max.mapping.bias |
        posteriors$start != posteriors$end

    posteriors$FLAGGED <- idx.ignore

    posteriors$log.ratio <- snv.lr[vcf.ids]
    posteriors$depth <- depth
    posteriors$prior.somatic <- prior.somatic[vcf.ids]
    posteriors$prior.contamination <- prior.cont[vcf.ids]
    posteriors$on.target <- info(vcf[vcf.ids])$OnTarget
    posteriors$seg.id <- queryHits(ov)

    if (!is.null(mapping.bias$pon.count)) {
        posteriors$pon.count <- mapping.bias$pon.count[vcf.ids]
        idx.ignore <- idx.ignore | 
            (posteriors$pon.count > max.pon & posteriors$prior.somatic > 0.1)
        posteriors$FLAGGED <- idx.ignore
    }
    # change seqnames to chr
    colnames(posteriors)[1] <- "chr"    

    ret <- list(
        llik = sum(log(rm.snv.posteriors[!idx.ignore])) - sum(idx.ignore), 
        likelihoods = likelihoods, 
        posteriors = posteriors, 
        vcf.ids = vcf.ids, 
        posterior.contamination = 0)

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
    # set sub-clonal (ML.M=0) to 1
    xx$ML.M[xx$ML.SOMATIC & !xx$ML.M] <- 1
    xx$POSTERIOR.SOMATIC <- apply(snv.posteriors[, grep("SOMATIC.M", 
        colnames(snv.posteriors))], 1, sum, na.rm=TRUE)
    xx[, c(1,3,2)]
}

.checkFraction <- function(x, name) {
    if (!is.numeric(x) || length(x) !=1 || 
        x < 0 || x > 1) {
        .stopUserError(name, " not within expected range or format.")
    }
}

.checkParameters <- function(test.purity, min.ploidy, max.ploidy, 
    max.non.clonal, max.homozygous.loss, sampleid, prior.K, 
    prior.contamination, prior.purity, iterations, min.gof, model.homozygous, 
    interval.file, log.ratio.calibration, test.num.copy) {
    if (min(test.purity) <= 0 || max(test.purity) > 0.99) 
        .stopUserError("test.purity not within expected range.")
    if (min.ploidy <= 0 || max.ploidy <= 2) 
        .stopUserError("min.ploidy or max.ploidy not within expected range.")

    if (min(test.num.copy) < 0) 
        .stopUserError("test.num.copy not within expected range.")
    
    if (min(test.num.copy) > 0 || max(test.num.copy)>8) 
        flog.warn("test.num.copy outside recommended range.")
                   
    .checkFraction(max.non.clonal, "max.non.clonal")
    .checkFraction(max.homozygous.loss[1], "max.homozygous.loss")
    .checkFraction(prior.K, "prior.K")
    .checkFraction(prior.contamination, "prior.contamination")
    .checkFraction(min.gof, "min.gof")

    tmp <- sapply(prior.purity, .checkFraction, "prior.purity")

    if (!is.null(sampleid) && (!is(sampleid, "character") ||
        length(sampleid) != 1)) {
        .stopUserError("sampleid not a character string.")
    }
    if (abs(1-sum(prior.purity)) > 0.02) {
        .stopUserError("prior.purity must add to 1. Sum is ", sum(prior.purity))
    }    
    if (length(prior.purity) != length(test.purity)) {
        .stopUserError("prior.purity must have the same length as ",
            "test.purity.")
    }    
    if (!is.null(interval.file) && !file.exists(interval.file)) {
        .stopUserError("interval.file ", interval.file, " not found.")
    }

    stopifnot(is.numeric(min.ploidy))
    stopifnot(is.numeric(max.ploidy))
    stopifnot(is.numeric(test.purity))
    stopifnot(is.numeric(iterations))
    stopifnot(is.numeric(log.ratio.calibration))
    stopifnot(is.logical(model.homozygous))

    if (iterations < 10 || iterations > 250) {
        .stopUserError("Iterations not in the expected range from 10 to 250.")
    }    
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
.filterDuplicatedResults <- function(results, purity.cutoff = 0.1) {
    if (length(results) < 2) 
        return(results)
    idx.duplicated <- rep(FALSE, length(results))

    for (i in seq_len(length(results)-1)) {
        for (j in seq(i+1,length(results))) {
            if (abs(results[[i]]$purity - results[[j]]$purity) < purity.cutoff &&
                abs(results[[i]]$ploidy - results[[j]]$ploidy) /
                results[[i]]$ploidy < 0.1) {
                idx.duplicated[j] <- TRUE
            } 
        }
    }
    results[!idx.duplicated]
}
.filterDuplicatedCandidates <- function(candidates) {
    if (nrow(candidates) < 2) 
        return(candidates)
    idx.duplicated <- rep(FALSE, nrow(candidates))

    for (i in seq_len(nrow(candidates)-1)) {
        for (j in seq(i+1,nrow(candidates))) {
            if (abs(candidates$purity[i] - candidates$purity[j]) < 0.1 &&
                abs(candidates$tumor.ploidy[i] - candidates$tumor.ploidy[j]) / 
                candidates$tumor.ploidy[i] < 0.1) {
                idx.duplicated[j] <- TRUE
            } 
        }
    }
    candidates[!idx.duplicated, ]
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
.isRareKaryotype <- function(ploidy) {
    ploidy > 4.5 || ploidy < 1.5
}
    
.flagResult <- function(result, max.non.clonal = 0.2, min.gof, 
    use.somatic.status, model.homozygous) {
    result$flag_comment <- NA
    result$flag <- .failedNonAberrant(result)
    if (result$flag) {
        result$flag_comment <- .appendComment(result$flag_comment, 
            "NON-ABERRANT")
    }
    if (result$fraction.subclonal > max.non.clonal*0.75) {
        result$flag <- TRUE
        result$flag_comment <- .appendComment(result$flag_comment, 
            "POLYGENOMIC")
    }
    if (result$fraction.homozygous.loss > 0.01) {
        result$flag <- TRUE
        result$flag_comment <- .appendComment(result$flag_comment, 
            "EXCESSIVE LOSSES")
    }
    if (.isRareKaryotype(result$ploidy)) {
        result$flag <- TRUE
        result$flag_comment <- .appendComment(result$flag_comment, 
            "RARE KARYOTYPE")
    }
    if (result$purity < 0.3) {
        result$flag <- TRUE
        result$flag_comment <- .appendComment(result$flag_comment, 
            "LOW PURITY")
    }
    if (result$purity > 0.9 && !model.homozygous && 
        (!is.null(use.somatic.status) && !use.somatic.status)) {
        result$flag <- TRUE
        result$flag_comment <- .appendComment(result$flag_comment, 
            "HIGH PURITY AND model.homozygous=FALSE")
    }
    result$GoF <- .getGoF(result)

    if (!is.null(result$GoF) && result$GoF < min.gof) {
        result$flag <- TRUE
        result$flag_comment <- .appendComment(result$flag_comment, 
            paste0("POOR GOF (", round(result$GoF*100,digits=1),"%)"))
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

.getFractionLoh <- function(result) {
    if (is.null(result$SNV.posterior)) return(0)
    pp <- result$SNV.posterior$posteriors
    x1 <- unique(pp$seg.id[pp$ML.M.SEGMENT==0])
    sum(result$seg$size[x1])/sum(result$seg$size)
}
.getGoF <- function(result) {
    if (is.null(result$SNV.posterior)) return(0)
    r <- result$SNV.posterior$posteriors
    e <- (r$ML.AR-r$AR.ADJUSTED)^2
    maxDist <- 0.2^2
    r2 <- max(1-mean(e,na.rm=TRUE)/maxDist,0)
    return(r2)
}    

.flagResults <- function(results, max.non.clonal = 0.2, max.logr.sdev, 
    logr.sdev, max.segments, min.gof, flag = NA, flag_comment = NA, 
    dropout=FALSE, use.somatic.status=TRUE, model.homozygous=FALSE) {
    if (length(results) < 1) return(results)

    results <- lapply(results, .flagResult, max.non.clonal=max.non.clonal, 
        min.gof=min.gof, use.somatic.status=use.somatic.status, 
        model.homozygous=model.homozygous)

    number.segments <- nrow(results[[1]]$seg)
    
    # some global flags
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
    
    if (!is.na(flag) && flag) {
        for (i in seq_along(results)) {
            results[[i]]$flag <- TRUE
            results[[i]]$flag_comment <- .appendComment(results[[i]]$flag_comment, 
                flag_comment)
        }
    }
    results
}
.getGeneCalls <- function(seg.adjusted, tumor, log.ratio, fun.focal, 
    args.focal, chr.hash) {
    args.focal <- c(list(seg = seg.adjusted), args.focal)
    focal <- do.call(fun.focal, args.focal)
    abs.gc <- GRanges(seqnames = .add.chr.name(seg.adjusted$chrom, chr.hash), IRanges(start = seg.adjusted$loc.start, 
        end = seg.adjusted$loc.end))

    # that will otherwise mess up the log-ratio means, mins and maxs
    idx <- which(!is.nan(log.ratio) & is.finite(log.ratio) & tumor$Gene != ".")
    if (!length(idx)) return(NA)
    tumor <- tumor[idx]
    log.ratio <- log.ratio[idx]
    
    ov <- findOverlaps(tumor, abs.gc)
    # use funky data.table to calculate means etc. in two lines of code.
    dt <- data.table(as.data.frame(tumor[queryHits(ov)]), C = seg.adjusted[subjectHits(ov), "C"], 
        seg.mean = seg.adjusted[subjectHits(ov), "seg.mean"], LR = log.ratio[queryHits(ov)], 
        seg.id = subjectHits(ov), seg.length = seg.adjusted$size[subjectHits(ov)], 
        focal = focal[subjectHits(ov)]
        )
    # some targets have multipe genes assigned?
    if (sum(grepl(",", dt$Gene))) {
        dt <- dt[, list(Gene = unlist(strsplit(as.character(Gene), ",", fixed=TRUE))), 
            by = list(seqnames, start, end, C, seg.mean, seg.id, seg.length, LR, focal)]
    }

    gene.calls <- data.frame(dt[, list(chr = seqnames[1], start = min(start), 
        end = max(end), C = as.double(median(C)), seg.mean = median(seg.mean),
        seg.id = seg.id[which.min(seg.length)], 
        .min.segid=min(seg.id), .max.segid=max(seg.id),
        number.targets = length(start), gene.mean = mean(LR, na.rm = TRUE), 
        gene.min = min(LR, na.rm = TRUE), gene.max = max(LR, na.rm = TRUE), 
        focal = focal[which.min(seg.length)]), by = Gene], row.names = 1)
    breakpoints <- gene.calls$.max.segid - gene.calls$.min.segid
    gene.calls$breakpoints <- breakpoints
    gene.calls    
}


.extractCountMatrix <- function(coverages) {
    useCounts <- .coverageHasCounts(coverages)
    if (useCounts) {
        return(do.call(cbind, 
            lapply(coverages, function(x) x$counts)))
    }
    # TODO: remove spring 2018 release
    flog.info("Coverage file does not contain read count information, %s", 
        "using total coverage for calculating log-ratios.")
    do.call(cbind, 
      lapply(coverages, function(x) x$coverage))
}
.coverageHasCounts <- function(coverages) {
    for (i in seq_along(coverages)) 
        if (sum(!is.na(coverages[[i]]$counts))==0) return(FALSE)
    return(TRUE)        
}
.checkSymbolsChromosome <- function(tumor, blank=c("", ".")) {
    if (is.null(tumor$Gene)) {
        tumor$Gene <- "."
        return(tumor)
    }    
    chrsPerSymbol <- lapply(split(seqnames(tumor), tumor$Gene), unique)
    nonUniqueSymbols <- names(chrsPerSymbol[sapply(chrsPerSymbol, length)>1])
    idx <- tumor$Gene %in% nonUniqueSymbols
    idxBlank <- tumor$Gene %in% blank
    tumor$Gene <- as.character(tumor$Gene)
    tumor$Gene[idx] <- paste(tumor$Gene, tumor$chr, sep=".")[idx]
    tumor$Gene[idxBlank] <- "."
    tumor
}
    
.get2DPurityGrid <- function(test.purity, by=1/30) {
    startPurity <- max(0.1, min(test.purity))
    endPurity <- min(0.99, max(test.purity))
    grid <- seq(startPurity, endPurity, by=by)   
    if (startPurity < 0.34 && endPurity > 0.35) {
        grid <- c(seq(startPurity, 0.34, by=1/50), 
        seq(0.35, endPurity, by=by))
    } 
    grid       
}
    
.optimizeGrid <- function(test.purity, min.ploidy, max.ploidy, test.num.copy = 0:7, 
    exon.lrs, seg, sd.seg, li, max.exon.ratio, max.non.clonal, BPPARAM) {
    ploidy.grid <- seq(min.ploidy, max.ploidy, by = 0.2)
    if (min.ploidy < 1.8 && max.ploidy > 2.2) {
        ploidy.grid <- c(seq(min.ploidy, 1.8, by = 0.2), 1.9, 2, 2.1, seq(2.2, max.ploidy, 
            by = 0.2))
    }

    .optimizeGridPurity <- function(p) {
        b <- 2 * (1 - p)
        log.ratio.offset <- 0
        lapply(ploidy.grid, function(D) {
            dt <- p/D
            llik.all <- lapply(seq_along(exon.lrs), function(i) .calcLlikSegmentExonLrs(exon.lrs[[i]], 
                log.ratio.offset, max.exon.ratio, sd.seg, dt, b, D, test.num.copy))
            subclonal <- vapply(llik.all, which.max, double(1)) == 1
            subclonal.f <- length(unlist(exon.lrs[subclonal]))/length(unlist(exon.lrs))
            if (subclonal.f > max.non.clonal) return(-Inf)
            sum(vapply(llik.all, max, double(1)))
        })
    }
    if (is.null(BPPARAM)) {
        mm <- lapply(test.purity, .optimizeGridPurity)
    } else {
        mm <- BiocParallel::bplapply(test.purity, .optimizeGridPurity, BPPARAM = BPPARAM)
    }
    mm <- sapply(mm, function(x) unlist(x))
    colnames(mm) <- test.purity
    rownames(mm) <- ploidy.grid

    if (!sum(as.vector(is.finite(mm)))) {
        .stopUserError("Cannot find valid purity/ploidy solution. ", 
            "This happens when input segmentations are garbage, most likely ",
            "due to a catastrophic sample QC failure. Re-check standard QC ",
            "metrics for this sample.")
    }

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
        
        if (sum(candidates$purity>grid[i] & 
            candidates$purity< grid[i+1], na.rm=TRUE) < 1) next

        # Nothing close to diplod in this range? Then add.
        if (min(abs(2 - candidates$tumor.ploidy[candidates$purity>grid[i] & 
            candidates$purity< grid[i+1]])) > 0.3) {

            mm.05 <- mm[, seq(t1+1, t2), drop=FALSE]
            
            # Find row most similar to normal diploid
            diploidRowId <- which.min(abs(2-as.numeric(row.names(mm.05))))
            # assert that rownames are still what they should be
            if (diploidRowId != which.min(abs(2-ploidy.grid))) {
                .stopRuntimeError("Cannot find diploid row in grid search.")
            }
            candidates <- rbind(candidates, 
                c(2, as.numeric(names(which.max(mm.05[diploidRowId, ]))), 
                max(mm.05[diploidRowId, ]), 2))

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
    candidates <- .filterDuplicatedCandidates(candidates)
    
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
    if (length(results) < 2) return(0)
    cs <- sapply((0:7)[-3], function(i) sapply(results, function(y)
                    sum(y$seg$size[y$seg$C == i])/sum(y$seg$size)))
    complexity <- apply(cs,1, function(z) sum(z>0.001))
    n <- sum(results[[1]]$seg$num.mark, na.rm=TRUE)
    -complexity*log(n)
}

.rankResults <- function(results) {
    if (length(results) < 1) return(results)  
    complexity <- .calcComplexityCopyNumber(results) 
    for (i in seq_along(results)) {
        if (is.null(results[[i]]$SNV.posterior)) {
            results[[i]]$total.log.likelihood <- results[[i]]$log.likelihood
        } else {
            results[[i]]$total.log.likelihood <- results[[i]]$log.likelihood/2 + results[[i]]$SNV.posterior$llik
        }
        results[[i]]$total.log.likelihood <- results[[i]]$total.log.likelihood + complexity[i]
    }
    idx.opt <- order(sapply(results, function(z) z$total.log.likelihood), decreasing = TRUE)
    
    results <- results[idx.opt]
    
    # remove solutions with -inf likelihood score
    results[!is.infinite(sapply(results, function(z) z$total.log.likelihood))]
}

.sampleOffsetFast <- function(test.num.copy, seg, exon.lrs, sd.seg, p, C, total.ploidy, 
    max.exon.ratio, simulated.annealing, log.ratio.calibration) {
    # Gibbs sample offset
    test.offset <- seq(sd.seg * -log.ratio.calibration, sd.seg * log.ratio.calibration, 
        by = 0.01)
    test.offset <- seq(p * -log.ratio.calibration, p * log.ratio.calibration * 0.2, length.out=12)

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
.sampleError <- function(subclonal, seg, exon.lrs, sd.seg, p, C, total.ploidy, max.exon.ratio, 
    simulated.annealing, iter, log.ratio.calibration, log.ratio.offset) {
    # Gibbs sample error
    test.error <- seq(sd.seg, sd.seg + sd.seg * 0.2, length.out = 5)

    seg.ids <- seq_len(nrow(seg))
    
    px.error <- lapply(test.error, function(error) vapply(seg.ids, 
        function(i) .calcLlikSegment(subclonal[i], exon.lrs[[i]] + log.ratio.offset[i], error, 
            p, C[i], total.ploidy, max.exon.ratio), double(1))
    )
    px.error.s <- sapply(px.error, sum, na.rm = TRUE)
    if (simulated.annealing) 
        px.error.s <- px.error.s * exp(iter/4)
    px.error.s <- exp(px.error.s - max(px.error.s))
    error <- test.error[min(which(runif(n = 1, min = 0, max = sum(px.error.s)) <= 
        cumsum(px.error.s)))]
}

.sampleOffset <- function(subclonal, seg, exon.lrs, sd.seg, p, C, total.ploidy, max.exon.ratio, 
    simulated.annealing, iter, log.ratio.calibration = 0.25) {
    # Gibbs sample offset
    test.offset <- seq(sd.seg * -log.ratio.calibration, sd.seg * log.ratio.calibration, 
        by = 0.01)
    test.offset <- seq(p * -log.ratio.calibration, p * log.ratio.calibration * 0.2, length.out=12)
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
.calcPuritySomaticVariants <- function(vcf, prior.somatic, tumor.id.in.vcf) {
    median(unlist(geno(vcf[prior.somatic > 0.5])$FA[, tumor.id.in.vcf]), na.rm = TRUE)/0.48
}
.createFakeLogRatios <- function(tumor, seg.file, sampleid, chr.hash, 
    model.homozygous=FALSE, max.logr.sdev) {
    if (!is.null(tumor$log.ratio)) {
         # calculate log.ratio sd in chunks of size 25 to estimate the 
         # segmented sd
         .robustSd <- function(d, size = 25) median(
            sapply(split(d, ceiling(seq_along(d) / size)), sd, na.rm = TRUE), 
            na.rm = TRUE)

        if (.robustSd(tumor$log.ratio) < max.logr.sdev) {
            flog.info("Found log2-ratio in tumor coverage data.")
            return(tumor$log.ratio)
        } else {
            flog.info("Provided log2-ratio looks too noisy, using segmentation only.")
        }
    }
    if (is(seg.file, "character")) {
        seg <- readSegmentationFile(seg.file, sampleid, model.homozygous = model.homozygous, verbose=FALSE)    
    } else {
        seg <-.checkSeg(seg.file, sampleid, model.homozygous, verbose=FALSE)
    }    
    seg.gr <- GRanges(seqnames = .add.chr.name(seg$chrom, chr.hash), 
                IRanges(start = round(seg$loc.start), end = seg$loc.end))

    ov <- findOverlaps(tumor, seg.gr)
    log.ratio <- seg$seg.mean[subjectHits(ov)]
    # sanity check, so that every exon has exactly one segment log-ratio
    log.ratio <- log.ratio[match(seq_len(length(tumor)), queryHits(ov))]
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
    chr.hash <- genomeStyles(species="Homo_sapiens")
    chr.hash <- chr.hash[!chr.hash$circular,]
    id <- which.max(apply(chr.hash[,-(1:3)],2,function(x) sum(ls %in%x)))+3
    chr.hash <- data.frame(chr=chr.hash[,id], number=seq_len(nrow(chr.hash)),
        row.names=chr.hash[,id])
    if (sum(!ls %in% chr.hash[,1]) == 0) return(chr.hash)
    data.frame(chr=as.factor(ls), number=seq_along(ls), row.names=ls)
}
.stopUserError <- function(...) {
    msg <- paste(c(...), collapse="")
    msg <- paste0(msg, "\n\nThis is most likely a user error due to",
        " invalid input data or parameters (PureCN ", 
        packageVersion("PureCN"), ").")
    flog.fatal(paste(strwrap(msg),"\n"))
    stop(paste(strwrap(msg),"\n"), call.= FALSE)
}
.stopRuntimeError <- function(...) {
    msg <- paste(c(...), collapse="")
    msg <- paste0(msg, "\n\nThis runtime error might be caused by",
        " invalid input data or parameters. Please report bug (PureCN ", 
        packageVersion("PureCN"), ").")
    flog.fatal(paste(strwrap(msg),"\n"))
    stop(paste(strwrap(msg),"\n"), call.= FALSE)
}
.logHeader <- function(l) {
    flog.info(strrep("-", 60))
    flog.info("PureCN %s", as.character(packageVersion("PureCN")))
    flog.info(strrep("-", 60))
    # which arguments are printable in a single line?
    idxSupported <- sapply(l, function(x) class(eval.parent(x))) %in% 
        c("character", "numeric", "NULL", "list", "logical") &
        sapply(l, function(x) object.size(eval.parent(x))) < 1024
    idxSmall <-     
        sapply(l[idxSupported], function(x) length(unlist(eval.parent(x)))) < 12
    idxSupported[idxSupported] <- idxSmall    

    l <- c(l[idxSupported],lapply(l[!idxSupported], function(x) "<data>"))
    argsStrings <- paste(sapply(seq_along(l), function(i) paste0("-", 
        names(l)[i], " ", paste(eval.parent(l[[i]]),collapse=","))),
        collapse=" ")
    flog.info("Arguments: %s", argsStrings)
}
.logFooter <- function() {
    flog.info("Done.")
    flog.info(strrep("-", 60))
}    
.calcGCmetric <- function(gc_bias, coverage, on.target) {
    idx <- which(coverage$on.target==on.target)
    if (!length(idx)) return(NA)
    gcbins <- split(coverage[idx]$average.coverage, 
        gc_bias[idx] < 0.5)
    mean(gcbins[[1]], na.rm=TRUE) / mean(gcbins[[2]], na.rm=TRUE) 
}
.checkGCBias <- function(normal, tumor, max.dropout, on.target=TRUE) {

    gcMetricNormal <- .calcGCmetric(tumor$gc_bias, normal, on.target)
    gcMetricTumor <- .calcGCmetric(tumor$gc_bias, tumor, on.target)

    if (is.na(gcMetricTumor)) return(FALSE)

    flog.info("AT/GC dropout%s: %.2f (tumor), %.2f (normal). ", 
        ifelse(on.target,""," (off-target regions)"), gcMetricTumor, gcMetricNormal)
    if (gcMetricNormal < max.dropout[1] || 
        gcMetricNormal > max.dropout[2] ||
        gcMetricTumor  < max.dropout[1] ||
        gcMetricTumor  > max.dropout[2]) {
        flog.warn("High GC-bias in normal or tumor. Is data GC-normalized?")
        return(TRUE)
    }
    return(FALSE)
}

.gcGeneToCoverage <- function(interval.file, min.coverage) {
    gc.data <- readCoverageFile(interval.file)
    gc.data$average.coverage <- min.coverage
    gc.data$coverage <- min.coverage * width(gc.data)
    gc.data
}

.getCentromerePositions <- function(centromeres, genome, style=NULL) {
    if (is.null(centromeres)) {
        data(centromeres, envir = environment())
        if (genome %in% names(centromeres)) {
            centromeres <- centromeres[[genome]]
            if (!is.null(style)) seqlevelsStyle(centromeres) <- style[1]
        } else {
            centromeres <- NULL
        }
    }
    centromeres
}
.checkArgs <- function(args, fname) {
    dups <- duplicated(names(args)) 
    if (sum(dups)) {
        args <- args[!dups]
        flog.warn("Duplicated arguments in %s", fname)
    }
    args
}
.calcFractionBalanced <- function(p) {
    sum(p$ML.C - p$ML.M.SEGMENT == p$ML.M.SEGMENT, na.rm=TRUE)/nrow(p)
}

# function to adjust log-ratios to segment mean
.postprocessLogRatios <- function(exon.lrs, seg.mean) {

    exon.lrs <- lapply(seq_along(exon.lrs), function(i) {
        if (length(exon.lrs[[i]]) > 3) return(exon.lrs[[i]])
        exon.lrs[[i]] - mean(exon.lrs[[i]]) + seg.mean[i]
    })

    return(exon.lrs)
}

.estimateContamination <- function(pp,  max.mapping.bias=NULL, min.fraction.chromosomes=0.8) {
    if (is.null(max.mapping.bias)) max.mapping.bias=0

    idx <- pp$GERMLINE.CONTHIGH + pp$GERMLINE.CONTLOW > 0.5 & 
        pp$MAPPING.BIAS >= max.mapping.bias

    if (!length(which(idx))) return(0)
    df <- data.frame(
        chr=pp$chr[idx], 
        AR=sapply(pp$AR.ADJUSTED[idx], function(x) ifelse(x>0.5, 1-x,x)),
        HIGHLOW=ifelse(pp$GERMLINE.CONTHIGH>pp$GERMLINE.CONTLOW, 
            "HIGH", "LOW")[idx]
    )
    # take the chromosome median and then average. the low count
    # might be biased in case contamination rate is < AR cutoff
    estimatedRate <- weighted.mean( 
        sapply(split(df$AR, df$HIGHLOW), median), 
        sapply(split(df$AR, df$HIGHLOW), length)
    )
    fractionChrs <- sum(unique(pp$chr) %in% df$chr)/length(unique(pp$chr))
    estimatedRate <- if (fractionChrs >= min.fraction.chromosomes) estimatedRate else 0
    estimatedRate
}

.calculate_ccf <- function(vaf, depth, purity, C){
    # see DOI: 10.1126/scitranslmed.aaa1408
    if (is.na(vaf)) return(c(NA, NA, NA))
    possible_ccfs <- seq(0.01, 1, 0.01)
    possible_vafs <- (purity * possible_ccfs)/
        ((2 * (1 - purity)) + (purity * C)) #Expected VAF for each CCF
    possible_vafs <- pmax(pmin(possible_vafs, 1), 0)    
    probs <- dbinom(x=round(vaf*depth), size = depth, prob = possible_vafs) #Prob of observed VAF
    names(probs) <- possible_ccfs
    if (!sum(probs)) {
        if (vaf > max(possible_vafs)) return(c(1, 1, 1))
        return(c(NA, NA, NA))
    }
    probs_norm <- probs / sum(probs) #Normalise to get posterior distribution
    probs_sort <- sort(probs_norm, decreasing = TRUE)
    probs_cum <- cumsum(probs_sort)
    n <- sum(probs_cum < 0.95) + 1 #Get 95% confidence interval (95% of probability)
    threshold <- probs_sort[n]
    cint  <- probs[probs_norm >= threshold]
    ccf_point <- as.numeric(names(which.max(probs_norm)))
    ccf_lower <- as.numeric(names(cint)[1])
    ccf_upper <- as.numeric(names(cint)[length(cint)])
    return(c(ccf_point, ccf_lower, ccf_upper))
}

.getExonLrs <- function(seg.gr, tumor, log.ratio, idx = NULL) {
    if (!is.null(idx)) {
        tumor <- tumor[idx]
        log.ratio <- log.ratio[idx]
    }    
    ov.se <- findOverlaps(seg.gr, tumor)
    exon.lrs <- lapply(seq_len(length(seg.gr)), function(i) log.ratio[subjectHits(ov.se)[queryHits(ov.se) == 
        i]])
    exon.lrs <- lapply(exon.lrs, function(x) subset(x, !is.na(x) & !is.infinite(x)))
    # estimate stand. dev. for target logR within targets. this will be used as proxy
    # for sample error.
    targetsPerSegment <- vapply(exon.lrs, length, integer(1))
    if (!sum(targetsPerSegment > 50, na.rm = TRUE) && !is.null(idx)) {
        .stopRuntimeError("Only tiny segments.")
    }
    exon.lrs
}
