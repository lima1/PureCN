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
.calcMsSegmentC <- function(yy, test.num.copy, Ci, prior.K) {
    max.M <- floor(Ci/2)
    idx.germline <- test.num.copy+length(test.num.copy)+1
    idx.somatic <- test.num.copy+1
    yys <- lapply(0:max.M, function(Mi) {
        for (i in test.num.copy) {
            n.cases.germ <- ifelse(Mi==Ci-Mi,1,2)
            tmp <- unique(c(0, 1,Mi, Ci-Mi))
            #tmp <- tmp[tmp>0]
            n.cases.somatic <- length(tmp)
            Cx <- max(i,Ci)

            if (i!=Mi && i!=Ci - Mi) {
                yy[,idx.germline[i+1]] <- yy[,idx.germline[i+1]] + log(1-prior.K) - log(Cx + 1 - n.cases.germ)

                # allow somatic mutations always have M=1
                if (i<=1) {
                    yy[,idx.somatic[i+1]] <- yy[,idx.somatic[i+1]] + log(prior.K) - log(n.cases.somatic)
                } else {
                    #message(paste(i, Mi, Ci, max.M, n.cases.somatic, n.cases.germ, Cx))
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
    best <- which.max(vapply(yys, function(x) sum(apply(x, 1, max)),double(1)))
    list(best=yys[[best]], M=test.num.copy[best])
}

# calculate likelihood scores for all possible minor copy numbers
.calcMsSegment <- function(xxi, test.num.copy, opt.C, prior.K) {
    lapply(seq_along(xxi), function(i).calcMsSegmentC( xxi[[i]], test.num.copy,
c(test.num.copy, round(opt.C))[i], prior.K))
}    

.calcSNVLLik <- function(vcf, tumor.id.in.vcf, ov, p, test.num.copy, 
    C.likelihood, C, opt.C, median.C, snv.model, prior.somatic, mapping.bias, snv.lr, 
    sampleid = NULL, cont.rate = 0.01, prior.K, max.coverage.vcf, non.clonal.M,
    model.homozygous=FALSE, error=0.001, max.mapping.bias=0.8, max.pon) {
    
    .dbeta <- function(x, shape1, shape2, log, size) dbeta(x=x, 
        shape1=shape1, shape2=shape2, log=log)
    if (snv.model=="betabin") {
        .dbeta <- function(x, shape1, shape2, log, size) {
            size <- pmin(size, max.coverage.vcf)
            dbetabinom.ab(x=round(x*size),shape1=shape1, shape2=shape2, 
                size=size, log=TRUE)
         }   
    }    

    prior.cont <- ifelse(info(vcf)$DB, cont.rate, 0)
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
    sd.ar <- sd(unlist(geno(vcf)$FA[, tumor.id.in.vcf]))

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

        lapply(seq(ncol(C.likelihood)), function(k) {
            Ci <- c(test.num.copy, opt.C[i])[k]       
            priorM <- log(max(Ci,1) + 1 + haploid.penalty)
            
            skip <- test.num.copy > Ci | C.likelihood[i, k] <=0

            p.ar <- lapply(c(0, 1), function(g) {
                cns <- test.num.copy
                if (!g) cns[1] <- non.clonal.M
                dbetax <- (p * cns + g * (1-p))/(p * Ci + 2 * (1-p))
                l <- lapply(seq_along(test.num.copy), function(j) {
                    if (skip[j]) return(mInf_all)
                    .dbeta(x = dbetax[j], 
                    shape1 = shape1,  
                    shape2 = shape2, log = TRUE, size=dp_all[idx]) - priorM
                })
                do.call(cbind,l)
            })
            
            p.ar.cont.1 <- .dbeta(x = (p * Ci + 2 * (1 - p - cont.rate))/
                  (p * Ci + 2 * (1 - p)),shape1=shape1, shape2=shape2, 
                  log=TRUE, size=dp_all[idx]) - priorM

            p.ar.cont.2 <- .dbeta(x = (cont.rate)/(p * Ci + 2 * (1 - p)), 
                  shape1=shape1, shape2=shape2,log=TRUE,size=dp_all[idx]) -
                                  priorM

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
        })
        
    })
    
    tmp <- lapply(seq_along(xx),function(i) .calcMsSegment(xx[[i]], 
               test.num.copy, opt.C[seg.idx[i]], prior.K))

    xx <- lapply(tmp, lapply, function(x) x$best)
    
    # Get segment M's for each SNV
    segment.M <- sapply(seq_along(tmp), function(i) 
        tmp[[i]][[min(C[seg.idx[i]], max(test.num.copy))+1]]$M)
    segment.M <- unlist(sapply(seq_along(seg.idx), function(i) 
        rep(segment.M[i], sum(seg.idx[i]==queryHits(ov)))))

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
    
    # this just adds a lot of helpful info to the SNV posteriors
    posteriors <- likelihoods/rowSums(likelihoods)
    xx <- .extractMLSNVState(posteriors)
    
    posteriors <- cbind(
        as.data.frame(rowRanges(vcf[vcf.ids]))[, 1:3], 
        posteriors, 
        xx, 
        ML.C = C[queryHits(ov)],
        ML.M.SEGMENT=segment.M
    )
    
    posteriors$ML.AR <- (p * posteriors$ML.M + 
        ifelse(posteriors$ML.SOMATIC, 0, 1) * 
        (1 - p))/(p * posteriors$ML.C + 2 * (1 - p))
    posteriors$ML.AR[ posteriors$ML.AR > 1] <- 1 

    posteriors$AR <- unlist(geno(vcf[vcf.ids])$FA[, tumor.id.in.vcf])
    posteriors$AR.ADJUSTED <- posteriors$AR / mapping.bias$bias[vcf.ids]
    posteriors$AR.ADJUSTED[posteriors$AR.ADJUSTED>1] <- 1
    posteriors$MAPPING.BIAS <- mapping.bias$bias[vcf.ids]
    # Extract LOH
    posteriors$ML.LOH <- (posteriors$ML.M == posteriors$ML.C | 
        posteriors$ML.M == 0 | posteriors$ML.C == 1)
    
    posteriors$CN.SUBCLONAL <- subclonal
    posteriors$CELLFRACTION <-  (posteriors$AR/posteriors$ML.M)*
        (p*posteriors$ML.C+2*(1-p))/p
    posteriors$CELLFRACTION[!posteriors$ML.SOMATIC] <- NA

    rm.snv.posteriors <- apply(likelihoods, 1, max) 
    idx.ignore <- rm.snv.posteriors == 0 | 
        posteriors$MAPPING.BIAS < max.mapping.bias
    posteriors$FLAGGED <- idx.ignore    

    posteriors$log.ratio <- snv.lr[vcf.ids]
    posteriors$depth <-as.numeric(geno(vcf[vcf.ids])$DP[, tumor.id.in.vcf])
    posteriors$prior.somatic <- prior.somatic[vcf.ids]
    posteriors$prior.contamination <- prior.cont[vcf.ids]
    posteriors$on.target <- info(vcf[vcf.ids])$OnTarget
    posteriors$seg.id <- queryHits(ov)


    if (!is.null(mapping.bias$pon.count)) {
        posteriors$pon.count <- mapping.bias$pon.count[vcf.ids]
        idx.ignore <- idx.ignore | 
            ( posteriors$pon.count > max.pon & !info(vcf[vcf.ids])$DB )
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
    gc.gene.file) {
    if (min(test.purity) <= 0 || max(test.purity) > 0.99) 
        .stopUserError("test.purity not within expected range.")
    if (min.ploidy <= 0 || max.ploidy <= 2) 
        .stopUserError("min.ploidy or max.ploidy not within expected range.")

    .checkFraction(max.non.clonal, "max.non.clonal")
    .checkFraction(max.homozygous.loss[1], "max.homozygous.loss")
    .checkFraction(prior.K, "prior.K")
    .checkFraction(prior.contamination, "prior.contamination")
    .checkFraction(min.gof, "min.gof")

    tmp <- sapply(prior.purity, .checkFraction, "prior.purity")

    if (!is.null(sampleid) && ( class(sampleid) != "character" ||
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
    if (!is.null(gc.gene.file) && !file.exists(gc.gene.file)) {
        .stopUserError("gc.gene.file ", gc.gene.file, " not found.")
    }

    stopifnot(is.numeric(min.ploidy))
    stopifnot(is.numeric(max.ploidy))
    stopifnot(is.numeric(test.purity))
    stopifnot(is.numeric(iterations))
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
.filterDuplicatedCandidates <- function(candidates) {
    if (nrow(candidates) < 2) 
        return(candidates)
    idx.duplicated <- rep(FALSE, nrow(candidates))

    for (i in seq_len(nrow(candidates)-1)) {
        for (j in seq(i+1,nrow(candidates)) ) {
            if ( abs( candidates$purity[i] - candidates$purity[j] ) < 0.1 &&
                 abs( candidates$tumor.ploidy[i] - candidates$tumor.ploidy[j] ) / candidates$tumor.ploidy[i] < 0.1) {
                idx.duplicated[j] <- TRUE
            }    
        }    
    }
    candidates[!idx.duplicated,]
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
    if (result$ploidy > 4.5 || result$ploidy < 1.5) {
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
            paste0("POOR GOF (", round(result$GoF*100,digits=1),"%)") )
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
    idx <- which(!is.nan(log.ratio) & !is.infinite(log.ratio) & tumor$Gene != ".")
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

.addVoomToGeneCalls <- function(results, tumor.coverage.file, normalDB, 
    gc.gene.file) {
    y <- .voomTargets(tumor.coverage.file, normalDB, gc.gene.file,
        num.normals=10, plot.voom=FALSE)

    logFC <- y$coefficients[rownames(results[[1]]$gene.calls),2]
    p <- y$p.value[rownames(results[[1]]$gene.calls),2]

    for (i in seq_along(results)) {
        results[[i]]$gene.calls$.voom.gene.mean <- logFC
        results[[i]]$gene.calls$pvalue <- p
   }
   results
}

.voomTargets <- function(tumor.coverage.file, normalDB, gc.gene.file, 
    num.normals=NULL, plot.voom=FALSE) {
    
    countMatrix <- .voomCountMatrix(tumor.coverage.file, normalDB=normalDB, num.normals=num.normals)
    gc.data <- read.delim(gc.gene.file, as.is=TRUE)
    gc.pos <- .checkSymbolsChromosome(GRanges(gc.data[,1], Gene=gc.data$Gene))

    geneCountMatrix <- do.call(rbind, lapply(split(data.frame(countMatrix), gc.pos$Gene),
         apply, 2, sum, na.rm=TRUE))
    
    geneCountMatrix <- geneCountMatrix[complete.cases(geneCountMatrix),]
    dge <- edgeR::DGEList(geneCountMatrix, 
            group=c("tumor", rep("normal", ncol(countMatrix)-1)))
    v <- limma::voomWithQualityWeights(dge, plot=plot.voom)
    y <- limma::lmFit(v, design = stats::model.matrix(~dge$samples$group))
    y <- limma::eBayes(y)
}

.voomLogRatio <- function(tumor.coverage.file, normal.coverage.files=NULL, 
    normalDB=NULL, num.normals=NULL, plot.voom=TRUE) {
    countMatrix <- .voomCountMatrix(tumor.coverage.file, normal.coverage.files, 
                                    normalDB, num.normals)
    
    idx <- complete.cases(countMatrix)
    dge <- edgeR::DGEList(countMatrix[idx,], 
            group=c("tumor", rep("normal", ncol(countMatrix)-1)))
    v <- limma::voomWithQualityWeights(dge, plot=plot.voom)
    y <- limma::lmFit(v, design = stats::model.matrix(~dge$samples$group))
    y <- limma::eBayes(y)
    logRatio <- rep(NA, nrow(countMatrix))
    logRatio[idx] <- y$coefficients[,2]
    logRatio
}

.voomCountMatrix <- function(tumor.coverage.file, normal.coverage.files=NULL, 
                             normalDB, num.normals) {
    if (is.null(normal.coverage.files)) {
        if (!is.null(num.normals)) {
            normal.coverage.files <- findBestNormal(tumor.coverage.file, 
                normalDB, num.normals=num.normals)
        } else {
            normal.coverage.files <- normalDB[[1]]
        }    
    }
    if (is.character(tumor.coverage.file)) {
        tumor <- readCoverageFile(tumor.coverage.file)
    } else {
        tumor <- tumor.coverage.file
    }
    normals <- .readNormals(normal.coverage.files)

    countMatrix <- do.call(cbind, 
        lapply(c(list(tumor), normals), function(x) x$average.coverage))
}    

.checkSymbolsChromosome <- function(tumor, blank=c("", ".")) {
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
    exon.lrs, seg, sd.seg, li, max.exon.ratio, max.non.clonal) {
    ploidy.grid <- seq(min.ploidy, max.ploidy, by = 0.2)
    if (min.ploidy < 1.8 && max.ploidy > 2.2) {
        ploidy.grid <- c(seq(min.ploidy, 1.8, by = 0.2), 1.9, 2, 2.1, seq(2.2, max.ploidy, 
            by = 0.2))
    }
    mm <- lapply(test.purity, function(p) {
        b <- 2 * (1 - p)
        log.ratio.offset <- 0
        lapply(ploidy.grid, function(D) {
            dt <- p/D
            llik.all <- lapply(seq_along(exon.lrs), function(i) .calcLlikSegmentExonLrs(exon.lrs[[i]], 
                log.ratio.offset, max.exon.ratio, sd.seg, dt, b, D, test.num.copy))
            subclonal <- vapply(llik.all, which.max, double(1)) == 1
            subclonal.f <- length(unlist(exon.lrs[subclonal]))/length(unlist(exon.lrs))
            if (subclonal.f > max.non.clonal) 
                return(-Inf)
            llik <- sum(vapply(llik.all, max, double(1)))
        })
    })
    mm <- sapply(mm, function(x) unlist(x))
    colnames(mm) <- test.purity
    rownames(mm) <- ploidy.grid

    if (!sum(as.vector(is.finite(mm))) ) {
        .stopUserError("Cannot find valid purity/ploidy solution. ", 
            "This happens when input segmentations are garbage.")
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
        
        if ( sum(candidates$purity>grid[i] & 
            candidates$purity< grid[i+1], na.rm=TRUE) < 1) next

        # Nothing close to diplod in this range? Then add.
        if (min(abs(2 - candidates$tumor.ploidy[candidates$purity>grid[i] & 
            candidates$purity< grid[i+1] ])) > 0.3) {

            mm.05 <- mm[,seq(t1+1,t2),drop=FALSE]
            
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
    cs <- sapply((0:7)[-3], function(i) sapply(results, function(y)
                    sum(y$seg$size[y$seg$C == i])/sum(y$seg$size)))
    complexity <- apply(cs,1, function(z) sum(z>0.001))
    n <- sum(results[[1]]$seg$num.mark, na.rm=TRUE)
    penalty <- -complexity*log(n)
}

.rankResults <- function(results) {
    if (length(results) < 1) return(results)  
    # max.ll <- max(sapply(results, function(z) z$log.likelihood)) max.snv.ll <-
    # max(sapply(results, function(z) z$SNV.posterior$llik))
    complexity <- .calcComplexityCopyNumber(results) 
    for (i in seq_along(results)) {
        if (is.null(results[[i]]$SNV.posterior)) {
            results[[i]]$total.log.likelihood <- results[[i]]$log.likelihood
        } else {
            # results[[i]]$total.log.likelihood <- (results[[i]]$log.likelihood/max.ll) +
            # 2*(results[[i]]$SNV.posterior$llik/max.snv.ll)
            results[[i]]$total.log.likelihood <- results[[i]]$log.likelihood/2 + results[[i]]$SNV.posterior$llik
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
.calcPuritySomaticVariants <- function(vcf, prior.somatic, tumor.id.in.vcf) {
    median(unlist(geno(vcf[prior.somatic > 0.5])$FA[, tumor.id.in.vcf]), na.rm = TRUE)/0.48
}
.loadSegFile <- function(seg.file, sampleid, model.homozygous=FALSE, 
    verbose=TRUE) {
    if (is.null(seg.file)) return(NULL)
    seg <- read.delim(seg.file)
    .checkSeg(seg, sampleid, model.homozygous, verbose)
}
.checkSeg <- function(seg, sampleid, model.homozygous, verbose=TRUE) {
    
    required.colnames <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", 
        "seg.mean")
    required.colnames2 <- c("ID", "chromosome", "start", "end", "num_probes", 
        "mean")
    if (ncol(seg) > length(required.colnames)) {
        seg <- seg[,1:length(required.colnames)]
    }    
    if (identical(colnames(seg), required.colnames2)) {
        colnames(seg) <- required.colnames
    }    
    
    # The smallest possible log-ratio is about 8
    # for 0.99 purity and high ploidy.
    # remove artifacts with lower log-ratio
    if (!model.homozygous && min(seg$seg.mean, na.rm=TRUE) < -8) {
        nBefore <- nrow(seg)
        seg <- seg[which(seg$seg.mean >= -8 | seg$num.mark >= 4),]
        if (verbose) flog.warn("Removing %i short segments with log-ratio < -8.", 
            nBefore-nrow(seg))
    }    

    if (!identical(colnames(seg), required.colnames)) {
        .stopUserError(paste("Segmentation file expected with colnames", 
                paste(required.colnames, collapse = ", ")))
    }
    segs <- split(seg, seg$ID)
    matchedSeg <- match(make.names(sampleid), make.names(names(segs)))

    if (length(segs)==1) {
        if (!is.null(sampleid) && is.na(matchedSeg)) {
            flog.warn("Provided sampleid (%s) does not match %s found in %s",
                      sampleid, names(segs)[1], "segmentation.")
        }   
        matchedSeg <- 1 
    } else if (is.null(sampleid)) {
        .stopUserError("seg.file contains multiple samples and sampleid missing.")
    } else if (is.na(matchedSeg)) {
        .stopUserError("seg.file contains multiple samples and sampleid does not match any.")
    } else {
        seg <- segs[[matchedSeg]]
    }    
    seg
}
       
.createFakeLogRatios <- function(tumor, seg.file, sampleid, chr.hash, 
    model.homozygous=FALSE) {
    if (class(seg.file)=="character") {
        seg <- .loadSegFile(seg.file, sampleid, model.homozygous, verbose=FALSE)    
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
#.ffpeCleanLogRatio <- function(log.ratio, window=20) {
#   dlr <- c(0, diff(log.ratio))
#   start <- seq(1,length(dlr), by=window)
#   end <- seq(window,length(dlr), by=window)
#   sds <- sapply(seq_along(end), function(i) sd(dlr[start[i]:end[i]], na.rm=TRUE))
##   ids <- which(sds > quantile(sds, na.rm=TRUE, p=1-0.001))
#}

.stopUserError <- function(...) {
    msg <- paste(c(...), collapse="")
    msg <- paste0(msg, "\n\nThis is most likely a user error due to",
        " invalid input data or parameters (PureCN ", 
        packageVersion("PureCN"), ").")
    flog.fatal(paste(strwrap(msg),"\n"))
    stop( paste(strwrap(msg),"\n"), call.= FALSE)
}
.stopRuntimeError <- function(...) {
    msg <- paste(c(...), collapse="")
    msg <- paste0(msg, "\n\nThis runtime error might be caused by",
        " invalid input data or parameters. Please report bug (PureCN ", 
        packageVersion("PureCN"), ").")
    flog.fatal(paste(strwrap(msg),"\n"))
    stop( paste(strwrap(msg),"\n"), call.= FALSE)
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
.calcGCmetric <- function(gc_bias, coverage) { 
    gcbins <- split(coverage$average.coverage, gc_bias < 0.5); 
    mean(gcbins[[1]], na.rm=TRUE)/mean(gcbins[[2]], na.rm=TRUE) 
}
.checkGCBias <- function(normal, tumor, max.dropout) {
    gcMetricNormal <- .calcGCmetric(tumor$gc_bias, normal)
    gcMetricTumor <- .calcGCmetric(tumor$gc_bias, tumor)
    flog.info("AT/GC dropout: %.2f (tumor), %.2f (normal). ", 
        gcMetricTumor, gcMetricNormal)
    if (gcMetricNormal < max.dropout[1] || 
        gcMetricNormal > max.dropout[2] ||
        gcMetricTumor  < max.dropout[1] ||
        gcMetricTumor  > max.dropout[2]) {
        flog.warn("High GC-bias in normal or tumor. Is data GC-normalized?")
        return(TRUE)
    }
    return(FALSE)
}

.gcGeneToCoverage <- function(gc.gene.file, min.coverage) {
    gc.data <- readCoverageFile(gc.gene.file)
    gc.data$average.coverage <- min.coverage
    gc.data$coverage <- min.coverage * width(gc.data)
    gc.data
}

.getCentromerePositions <- function(centromeres, genome) {
    if (is.null(centromeres)) {
        data(centromeres, envir = environment())
        if (genome %in% names(centromeres)) {
            centromeres <- centromeres[[genome]]
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

.getMajorityStateTargets <- function(ret, id, tumor, n=1) {
    seg <- ret$results[[id]]$seg
    states <- sapply(seq(0,7), function(i) sum(seg$num.mark[which(round(seg$C)==i)]))
    majorityC <- head(order(states,decreasing=TRUE),n)-1
    chr.hash <- .getChrHash(seqlevels(tumor))
    majorityGr <- GRanges(seqnames=.add.chr.name(seg$chrom, chr.hash), 
        IRanges(start=seg$loc.start, end=seg$loc.end))
    majorityGr <- majorityGr[seg$C %in% majorityC]
    idx <- overlapsAny(tumor, majorityGr)
    percentMajority <- sum(idx)/length(idx)*100
    if (percentMajority < 50 && n < 3) {
        return(.getMajorityStateTargets(ret,id,tumor,n+1))
    }
    flog.info("Majority Copy Number State%s: %s (%.1f%% of targets)", 
        ifelse(length(majorityC)>1, "s",""), 
        paste(majorityC, collapse=", "), 
        percentMajority )
    idx
}

# function to adjust log-ratios to segment mean
.postprocessLogRatios <- function(exon.lrs, seg.mean) {

    exon.lrs <- lapply(seq_along(exon.lrs), function(i) {
        if (length(exon.lrs[[i]]) > 3) return(exon.lrs[[i]])
        exon.lrs[[i]] - mean(exon.lrs[[i]]) + seg.mean[i]
    })

    return(exon.lrs)
}    
