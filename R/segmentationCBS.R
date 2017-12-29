#' CBS segmentation
#' 
#' The default segmentation function. This function is called via the
#' \code{fun.segmentation} argument of \code{\link{runAbsoluteCN}}.  The
#' arguments are passed via \code{args.segmentation}.
#' 
#' 
#' @param normal Coverage data for normal sample.
#' @param tumor Coverage data for tumor sample.
#' @param log.ratio Copy number log-ratios, one for each target in the coverage
#' files.
#' @param seg If segmentation was provided by the user, this data structure
#' will contain this segmentation. Useful for minimal segmentation functions.
#' Otherwise PureCN will re-segment the data. This segmentation function
#' ignores this user provided segmentation.
#' @param plot.cnv Segmentation plots.
#' @param sampleid Sample id, used in output files.
#' @param target.weight.file Can be used to assign weights to targets.
#' @param alpha Alpha value for CBS, see documentation for the \code{segment}
#' function.
#' @param undo.SD \code{undo.SD} for CBS, see documentation of the
#' \code{segment} function. If NULL, try to find a sensible default.
#' @param vcf Optional \code{CollapsedVCF} object with germline allelic ratios.
#' @param tumor.id.in.vcf Id of tumor in case multiple samples are stored in
#' VCF.
#' @param normal.id.in.vcf Id of normal in in VCF. Currently not used.
#' @param max.segments If not \code{NULL}, try a higher \code{undo.SD}
#' parameter if number of segments exceeds the threshold.
#' @param prune.hclust.h Height in the \code{hclust} pruning step. Increasing
#' this value will merge segments more aggressively. If NULL, try to find a
#' sensible default.
#' @param prune.hclust.method Cluster method used in the \code{hclust} pruning
#' step. See documentation for the \code{hclust} function.
#' @param chr.hash Mapping of non-numerical chromsome names to numerical names
#' (e.g. chr1 to 1, chr2 to 2, etc.). If \code{NULL}, assume chromsomes are
#' properly ordered.
#' @param centromeres A \code{GRanges} object with centromere positions.
#' Currently not supported in this function.
#' @return \code{data.frame} containing the segmentation.
#' @author Markus Riester
#' @references Olshen, A. B., Venkatraman, E. S., Lucito, R., Wigler, M. 
#' (2004). Circular binary segmentation for the analysis of array-based DNA 
#' copy number data. Biostatistics 5: 557-572.
#'
#' Venkatraman, E. S., Olshen, A. B. (2007). A faster circular binary 
#' segmentation algorithm for the analysis of array CGH data. Bioinformatics 
#' 23: 657-63.
#'
#' @seealso \code{\link{runAbsoluteCN}}
#' @examples
#' 
#' normal.coverage.file <- system.file("extdata", "example_normal.txt", 
#'     package="PureCN")
#' tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
#'     package="PureCN")
#' vcf.file <- system.file("extdata", "example_vcf.vcf", 
#'     package="PureCN")
#' interval.file <- system.file("extdata", "example_intervals.txt", 
#'     package="PureCN")
#' 
#' # The max.candidate.solutions, max.ploidy and test.purity parameters are set to
#' # non-default values to speed-up this example.  This is not a good idea for real
#' # samples.
#' ret <-runAbsoluteCN(normal.coverage.file=normal.coverage.file, 
#'     tumor.coverage.file=tumor.coverage.file, vcf.file=vcf.file, genome="hg19", 
#'     sampleid="Sample1", interval.file=interval.file, 
#'     max.candidate.solutions=1, max.ploidy=4, test.purity=seq(0.3,0.7,by=0.05), 
#'     fun.segmentation=segmentationCBS, args.segmentation=list(alpha=0.001))
#' 
#' @export segmentationCBS
#' @importFrom stats t.test hclust cutree
segmentationCBS <- function(normal, tumor, log.ratio, seg, plot.cnv, 
    sampleid, target.weight.file = NULL, alpha = 0.005, undo.SD =
        NULL, vcf = NULL, tumor.id.in.vcf = 1, normal.id.in.vcf = NULL,
    max.segments = NULL, prune.hclust.h = NULL, prune.hclust.method = "ward.D",
    chr.hash = NULL, centromeres = NULL) {

    if (is.null(chr.hash)) chr.hash <- .getChrHash(seqlevels(tumor))
    
    target.weights <- NULL
    if (!is.null(target.weight.file)) {
        target.weights <- read.delim(target.weight.file, as.is=TRUE)
        target.weights <- target.weights[match(as.character(tumor), 
            target.weights[,1]),2]
        flog.info("Target weights found, will use weighted CBS.")
        #if (!is.null(normal$se)) {
        #    flog.info("Log-ratio standard errors found, will use them in weights.")
        #    target.weights <- (target.weights + (1 / normal$se^2))/2
        #}     
    }
    x <- .CNV.analyze2(normal, tumor, log.ratio=log.ratio, plot.cnv=plot.cnv, 
        sampleid=sampleid, alpha=alpha, 
        weights=target.weights, sdundo=undo.SD, max.segments=max.segments,
        chr.hash=chr.hash) 
    if (!is.null(vcf)) {
        x <- .pruneByVCF(x, vcf, tumor.id.in.vcf, chr.hash=chr.hash)
        x <- .findCNNLOH(x, vcf, tumor.id.in.vcf, alpha=alpha, 
            chr.hash=chr.hash)
        x <- .pruneByHclust(x, vcf, tumor.id.in.vcf, h=prune.hclust.h, 
            method=prune.hclust.method, chr.hash=chr.hash)
    }
    idx.enough.markers <- x$cna$output$num.mark > 1
    rownames(x$cna$output) <- NULL
    x$cna$output[idx.enough.markers,]
}

.findCNNLOH <- function(x, vcf, tumor.id.in.vcf, alpha=0.005, min.variants=7, 
iterations=2, chr.hash) {
    for (iter in seq_len(iterations)) {
        seg <- x$cna$output
        seg.gr <- GRanges(seqnames=.add.chr.name(seg$chrom, chr.hash), 
            IRanges(start=seg$loc.start, end=seg$loc.end))
        ov <- findOverlaps(seg.gr, vcf)
        ar <- sapply(geno(vcf)$FA[,tumor.id.in.vcf], function(x) x[1])
        ar.r <- ifelse(ar>0.5, 1-ar, ar)

        segs <- split(seg, seq_len(nrow(seg)))
        foundCNNLOH <- FALSE
        for (i in seq_len(nrow(seg))) {
            sar <- ar.r[subjectHits(ov)][queryHits(ov)==i]
            if (length(sar) < 2 * min.variants) next
            min.variants.x <- max(min.variants, length(sar)*0.05)    
            bp <- which.min(sapply(seq(min.variants.x, length(sar)-min.variants.x, by=1), 
                function(i) sum(c(sd(sar[seq_len(i)]),
                    sd(sar[seq(i+1,length(sar))])))
            ))
            bp <- bp + min.variants.x-1
            x1 <- sar[seq_len(bp)]
            x2 <- sar[seq(bp+1,length(sar))]
            tt <- t.test(x1, x2, exact=FALSE)
            if ((abs(mean(x1) - mean(x2)) > 0.05 && tt$p.value < alpha) || 
                (abs(mean(x1) - mean(x2)) > 0.025 && tt$p.value < alpha && 
                min(length(x1), length(x2)) > min.variants * 3)
                ) {
                segs[[i]] <- rbind(segs[[i]], segs[[i]])            
                bpPosition <-  start(vcf[subjectHits(ov)][queryHits(ov)==i])[bp]
                segs[[i]]$loc.end[1] <- bpPosition
                segs[[i]]$loc.start[2] <- bpPosition+1
                foundCNNLOH <- TRUE    
            }
        }
        if (foundCNNLOH) {
            x$cna$output <- .updateNumMark(do.call(rbind, segs), x)
        } else {
            # no need for trying again if no CNNLOH found.
            break
        }   
    }
    x
}

# After merging, we need to figure out how many targets are covering the new 
# segments
.updateNumMark <- function(seg, x) {
    segGR <- GRanges(seqnames=seg$chrom, IRanges(start=seg$loc.start, 
        end=seg$loc.end))
    probeGR <- GRanges(seqnames=as.character(x$cna$data$chrom), 
        IRanges(start=x$cna$data$maploc, end=x$cna$data$maploc))
    probeToSegs <- findOverlaps(probeGR, segGR, select="first")
    seg$num.mark <- sapply(seq_len(nrow(seg)), function(i) sum(probeToSegs==i, 
        na.rm=TRUE))
    seg
}
        
.pruneByHclust <- function(x, vcf, tumor.id.in.vcf, h=NULL, method="ward.D", 
    min.variants=5, chr.hash, iterations=2) {
    for (iter in seq_len(iterations)) {
        seg <- x$cna$output
        seg.gr <- GRanges(seqnames=.add.chr.name(seg$chrom, chr.hash), 
            IRanges(start=seg$loc.start, end=seg$loc.end))
        ov <- findOverlaps(seg.gr, vcf)

        ar <- sapply(geno(vcf)$FA[,tumor.id.in.vcf], function(x) x[1])
        ar.r <- ifelse(ar>0.5, 1-ar, ar)
        dp <- geno(vcf)$DP[, tumor.id.in.vcf]

        xx <- sapply(seq_len(nrow(seg)) , function(i) {
            weighted.mean(
                ar.r[subjectHits(ov)][queryHits(ov) == i], 
                w=sqrt(dp[subjectHits(ov)][queryHits(ov) == i]),
                na.rm=TRUE)
        })

        if (is.null(h)) {
            h <- .getPruneH(seg)
            flog.info("Setting prune.hclust.h parameter to %f.", h)
        }   

        numVariants <- sapply(seq_len(nrow(seg)), function(i) 
            sum(queryHits(ov) == i))
        dx <- cbind(seg$seg.mean,xx)
        hc <- hclust(dist(dx), method=method)
        seg.hc <- data.frame(id=1:nrow(dx), dx, num=numVariants, 
            cluster=cutree(hc,h=h))[hc$order,]

        # cluster only segments with at least n variants    
        seg.hc <- seg.hc[seg.hc$num>=min.variants,]
        clusters <- lapply(unique(seg.hc$cluster), function(i) 
            seg.hc$id[seg.hc$cluster==i])
        clusters <- clusters[sapply(clusters, length)>1]
        
        x$cna$output$cluster.id <- NA
        for (i in seq_along(clusters)) {
            x$cna$output$cluster.id[clusters[[i]]] <- i
            x$cna$output$seg.mean[clusters[[i]]] <- 
                weighted.mean(x$cna$output$seg.mean[clusters[[i]]], 
                x$cna$output$num.mark[clusters[[i]]])
        }
        # merge consecutive segments with the same cluster id 
        merged <- rep(FALSE, nrow(x$cna$output))
        for (i in 2:nrow(x$cna$output)) {
            if (is.na(x$cna$output$cluster.id[i-1]) || 
                is.na(x$cna$output$cluster.id[i]) || 
                x$cna$output$chrom[i-1] != x$cna$output$chrom[i]  ||
                merged[i-1] ||
                x$cna$output$cluster.id[i-1] != x$cna$output$cluster.id[i]) next
                merged[i] <- TRUE

                x$cna$output$num.mark[i-1] <- x$cna$output$num.mark[i]+x$cna$output$num.mark[i-1]
                x$cna$output$size[i-1] <- x$cna$output$size[i]+x$cna$output$size[i-1]
                x$cna$output$loc.end[i-1] <- x$cna$output$loc.end[i]
        }
        x$cna$output <- x$cna$output[!merged,]
    }
    x
}
    

# looks at breakpoints, and if p-value is higher than max.pval, merge unless 
# there is evidence based on germline SNPs
.pruneByVCF <- function(x, vcf, tumor.id.in.vcf, min.size=5, max.pval=0.00001,
    iterations=3, chr.hash, debug=FALSE) {
    seg <- try(segments.p(x$cna), silent=TRUE)
    if (class(seg) == "try-error") return(x)
    for (iter in seq_len(iterations)) {
        seg.gr <- GRanges(seqnames=.add.chr.name(seg$chrom, chr.hash), 
            IRanges(start=seg$loc.start, end=seg$loc.end))
        ov <- findOverlaps(seg.gr, vcf)
        ar <- sapply(geno(vcf)$FA[,tumor.id.in.vcf], function(x) x[1])
        ar.r <- ifelse(ar>0.5, 1-ar, ar)
        merged <- rep(FALSE, nrow(seg))
        for (i in 2:nrow(seg)) {
            # don't try to merge chromosomes or very significant breakpoints
            if (is.na(seg$pval[i-1]) || seg$pval[i-1]<max.pval) next
            # don't merge when we have no germline data for segments    
            if (!(i %in% queryHits(ov) && (i-1) %in% queryHits(ov))) next
            ar.i <- list(
                ar.r[subjectHits(ov)][queryHits(ov)==i],
                ar.r[subjectHits(ov)][queryHits(ov)==i-1])
            if (length(ar.i[[1]]) < min.size || length(ar.i[[2]]) < min.size) next
            if (merged[i-1]) next
            
            p.t <- t.test(ar.i[[1]], ar.i[[2]], exact=FALSE)$p.value
            if (p.t>0.2) {
                merged[i] <- TRUE
                x$cna$output$seg.mean[i-1] <- weighted.mean(
                    c(seg$seg.mean[i],seg$seg.mean[i-1]),
                    w=c(seg$num.mark[i],seg$num.mark[i-1]))

                x$cna$output$num.mark[i-1] <- seg$num.mark[i]+seg$num.mark[i-1]
                x$cna$output$loc.end[i-1] <- seg$loc.end[i]
                seg$pval[i-1] <- seg$pval[i]
            }
            if (debug) message(paste(i, "LR diff:", 
                abs(seg$seg.mean[i]-seg$seg.mean[i-1]), "Size: ", 
                seg$num.mark[i-1], "PV:", p.t, "PV bp:",seg$pval[i-1], 
                "Merged:", merged[i],"\n", sep=" "))
        }
        x$cna$output <- x$cna$output[!merged,]
        seg <- seg[!merged,]
    }
    x
}

.getSDundo <- function(log.ratio, d=0.1) {
    # did user provide segmentation? Then the sd of the log-ratio
    # is wrong and shouldn't be used 
    isFakeLogRatio <- sum(abs(diff(log.ratio))<0.0001)/
        length(log.ratio)>0.9
    if (isFakeLogRatio) return(0)
    # a crude way of estimating purity - in heigher purity samples, we
    # can undo a little bit more aggressively 
    q <- quantile(log.ratio,p=c(d, 1-d))
    q.diff <- abs(q[1] - q[2])
    if (q.diff < 0.5) return(0.5)
    if (q.diff < 1) return(0.75)
    if (q.diff < 1.5) return(1)
    return(1.25)
}

.getPruneH <- function(seg, d=0.05) {
    seg <- seg[seg$num.mark>=1,]
    log.ratio <- unlist(lapply(seq_len(nrow(seg)), function(i) 
        rep(seg$seg.mean[i], seg$num.mark[i])))
    q <- quantile(log.ratio,p=c(d, 1-d))
    q.diff <- abs(q[1] - q[2])
    if (q.diff < 0.5) return(0.1)
    if (q.diff < 1) return(0.15)
    if (q.diff < 1.5) return(0.2)
    return(0.25)
}    

.fixStartEndPos <- function(segment.CNA.obj, normal) {
    segment.CNA.obj$output$loc.start <- start(normal[segment.CNA.obj$segRows$startRow])
    segment.CNA.obj$output$loc.end <- end(normal[segment.CNA.obj$segRows$endRow])
    segment.CNA.obj
}

# ExomeCNV version without the x11() calls 
.CNV.analyze2 <-
function(normal, tumor, log.ratio=NULL, weights=NULL, sdundo=NULL,
undo.splits="sdundo", alpha, eta=0.05, nperm=10000, sampleid=NULL,
plot.cnv=TRUE, max.segments=NULL, chr.hash=chr.hash) {
    
    if (is.null(sdundo)) {
        sdundo <- .getSDundo(log.ratio)
    }   
     
    CNA.obj <- .getCNAobject(log.ratio, normal, chr.hash, sampleid)
    
    # default parameters, then speed-up the segmentation by using precomputed
    # simulations 
    if (nperm == 10000 & alpha == 0.005 & eta == 0.05) {
        flog.info("Loading pre-computed boundaries for DNAcopy...")
        data(purecn.DNAcopy.bdry, package = "PureCN", 
              envir = environment())
        sbdry <- get("purecn.DNAcopy.bdry", envir = environment())
    }
    else {
        max.ones <- floor(nperm * alpha) + 1
        sbdry <- getbdry(eta, nperm, max.ones)
    }

    try.again <- 0
    while (try.again < 2) {
        flog.info("Setting undo.SD parameter to %f.", sdundo)
        if (!is.null(weights)) {
            # MR: this shouldn't happen. In doubt, count them as median.
            weights[is.na(weights)] <- median(weights, na.rm=TRUE)
            segment.CNA.obj <- segment(CNA.obj, 
                undo.splits = undo.splits, undo.SD = sdundo, sbdry = sbdry,
                nperm = nperm, verbose = 0, alpha = alpha, weights = weights)
        } else {
            segment.CNA.obj <- segment(CNA.obj, 
                undo.splits = undo.splits, undo.SD = sdundo, sbdry = sbdry,
                nperm = nperm, verbose = 0, alpha = alpha)
        }
        if (sdundo <= 0 || is.null(max.segments) || 
            nrow(segment.CNA.obj$output) < max.segments) break
        sdundo <- sdundo * 1.5
        try.again <- try.again + 1
    }

    segment.CNA.obj <- .fixStartEndPos(segment.CNA.obj, normal)

    if (plot.cnv) {
        plot(segment.CNA.obj, plot.type="s")
        plot(segment.CNA.obj, plot.type="w")
    }

    return(list(cna=segment.CNA.obj, logR=log.ratio))
}

.getSegSizes <- function(seg) {
    round(seg$loc.end-seg$loc.start+1)
}

.getCNAobject <- function(log.ratio, normal, chr.hash, sampleid) {
    CNA(log.ratio, .strip.chr.name(seqnames(normal), chr.hash),
        floor( (start(normal) + end(normal)) / 2), 
        data.type = "logratio", sampleid = sampleid)
}
