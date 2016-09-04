segmentationCBS <-
structure(function(# CBS segmentation
### The default segmentation function. This function is called via the 
### \code{fun.segmentation} argument of \code{\link{runAbsoluteCN}}. 
### The arguments are passed via \code{args.segmentation}.
##seealso<< \code{\link{runAbsoluteCN}}
normal, 
### GATK coverage data for normal sample.
tumor,  
### GATK coverage data for tumor sample.
log.ratio, 
### Copy number log-ratios, one for each target in the coverage files.
plot.cnv, 
### Segmentation plots.
min.coverage, 
### Minimum coverage in both normal and tumor.
sampleid=sampleid,
### Sample id, used in output files.
exon.weight.file=NULL,
### Can be used to assign weights to targets. Deprecated, use 
### \code{target.weight.file} instead.
target.weight.file=NULL,
### Can be used to assign weights to targets. 
alpha=0.005,
### Alpha value for CBS, see documentation for the \code{segment} function.
undo.SD=NULL,
### \code{undo.SD} for CBS, see documentation of the \code{segment} function.
### If NULL, try to find a sensible default.
vcf=NULL,
### Optional \code{CollapsedVCF} object with germline allelic ratios.
tumor.id.in.vcf=1,
### Id of tumor in case multiple samples are stored in VCF.
normal.id.in.vcf=NULL,
### Id of normal in in VCF. Currently not used.
max.segments=NULL,
### If not \code{NULL}, try a higher \code{undo.SD} parameter 
### if number of segments exceeds the threshold.
prune.hclust.h=NULL,
### Height in the \code{hclust} pruning step. Increasing this value 
### will merge segments more aggressively. If NULL, try to find a sensible
### default.
prune.hclust.method="ward.D",
### Cluster method used in the \code{hclust} pruning step. See
### documentation for the \code{hclust} function.
chr.hash=NULL,
### Mapping of non-numerical chromsome names to numerical names
### (e.g. chr1 to 1, chr2 to 2, etc.). If \code{NULL}, assume chromsomes
### are properly ordered.
verbose=TRUE
### Verbose output.
) {
    if (is.null(chr.hash)) chr.hash <- .getChrHash(tumor$chr)
    
    if (!is.null(exon.weight.file)) {
        message("exon.weight.file is deprecated, use ",
            "target.weight.file instead.")
        target.weight.file <- exon.weight.file
    }

    target.weights <- NULL
    if (!is.null(target.weight.file)) {
        target.weights <- read.delim(target.weight.file, as.is=TRUE)
        target.weights <- target.weights[match(as.character(tumor[,1]), 
            target.weights[,1]),2]
        if (verbose) message("Target weights found, will use weighted CBS.")
    }
    x <- .CNV.analyze2(normal, tumor, logR=log.ratio, plot.cnv=plot.cnv, 
        min.coverage=min.coverage, sampleid=sampleid, alpha=alpha, 
        weights=target.weights, sdundo=undo.SD, max.segments=max.segments,
        chr.hash=chr.hash, verbose=verbose) 
    if (!is.null(vcf)) {
        x <- .pruneByVCF(x, vcf, tumor.id.in.vcf, chr.hash=chr.hash)
        x <- .findCNNLOH(x, vcf, tumor.id.in.vcf, alpha=alpha, 
            chr.hash=chr.hash)
        x <- .pruneByHclust(x, vcf, tumor.id.in.vcf, h=prune.hclust.h, 
            method=prune.hclust.method, chr.hash=chr.hash, verbose=verbose)
    }
    idx.enough.markers <- x$cna$output$num.mark > 1
    rownames(x$cna$output) <- NULL
    if (verbose) {
        print(x$cna$output[idx.enough.markers,])
    }
    x$cna$output[idx.enough.markers,]
### \code{data.frame} containing the segmentation.    
},ex=function() {
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")
vcf.file <- system.file("extdata", "example_vcf.vcf", 
    package="PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
    package="PureCN")

# The max.candidate.solutions, max.ploidy and test.purity parameters are set to
# non-default values to speed-up this example.  This is not a good idea for real
# samples.
ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
    gatk.tumor.file=gatk.tumor.file, vcf.file=vcf.file, genome="hg19", 
    sampleid='Sample1', gc.gene.file=gc.gene.file, 
    max.candidate.solutions=1, max.ploidy=4, test.purity=seq(0.3,0.7,by=0.05), 
    fun.segmentation=segmentationCBS, args.segmentation=list(alpha=0.001))
})    

.findCNNLOH <- function( x, vcf, tumor.id.in.vcf, alpha=0.005, min.variants=7, 
iterations=2, chr.hash ) {
    for (iter in seq_len(iterations)) {
        seg <- x$cna$output
        seg.gr <- GRanges(seqnames=.add.chr.name(seg$chrom, chr.hash), 
            IRanges(start=seg$loc.start, end=seg$loc.end))
        ov <- findOverlaps(seg.gr, vcf)
        ar <- sapply(geno(vcf)$FA[,tumor.id.in.vcf], function(x) x[1])
        ar.r <- ifelse(ar>0.5, 1-ar, ar)

        segs <- split(seg, seq_len(nrow(seg)) )
        foundCNNLOH <- FALSE
        for (i in seq_len(nrow(seg)) ) {
            sar <-ar.r[subjectHits(ov)][queryHits(ov)==i]
            if (length(sar) < 2 * min.variants) next
            bp <- which.min(sapply(seq(min.variants, length(sar)-min.variants, by=1), 
                function(i) sum(c( sd(sar[seq_len(i)]),
                    sd(sar[seq(i+1,length(sar))])))
            ))
            bp <- bp + min.variants-1
            x1 <- sar[seq_len(bp)]
            x2 <- sar[seq(bp+1,length(sar))]
            tt <- t.test(x1, x2, exact=FALSE)
            if ( abs(mean(x1) - mean(x2)) > 0.05
                && tt$p.value < alpha) {
                segs[[i]] <- rbind(segs[[i]], segs[[i]])            
                bpPosition <-  start(vcf[subjectHits(ov)][queryHits(ov)==i])[bp]
                segs[[i]]$loc.end[1] <- bpPosition
                segs[[i]]$loc.start[2] <- bpPosition+1
                #message("Iteration: ", iter, " Found CNNLOH on ", seg$chrom[i], 
                #    " Means: ", mean(x1), " ", mean(x2), 
                #    " Lengths: ", length(x1), " ", length(x2)) 
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
    min.variants=5, chr.hash, iterations=2, verbose=TRUE) {
    for (iter in seq_len(iterations)) {
    seg <- x$cna$output
    #message("HCLUST: ", iter, " Num segment LRs: ", length(table(x$cna$output$seg.mean)))
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
        if (verbose) message("Setting prune.hclust.h parameter to ", h)
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

    for (i in seq_along(clusters)) {
        x$cna$output$seg.mean[clusters[[i]]] <- 
            weighted.mean(x$cna$output$seg.mean[clusters[[i]]], 
            x$cna$output$num.mark[clusters[[i]]])
    }
    }
    #message("HCLUST Num segment LRs: ", length(table(x$cna$output$seg.mean)))
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

.getSDundo <- function(log.ratio) {
    q <- quantile(log.ratio,p=c(0.1, 0.9))
    q.diff <- abs(q[1] - q[2])
    if (q.diff < 1) return(0.5)
    if (q.diff < 1.5) return(0.75)
    return(1)
}    

.getPruneH <- function(seg) {
    seg <- seg[seg$num.mark>=1,]
    log.ratio <- unlist(lapply(seq_len(nrow(seg)), function(i) 
        rep(seg$seg.mean[i], seg$num.mark[i])))
    q <- quantile(log.ratio,p=c(0.1, 0.9))
    q.diff <- abs(q[1] - q[2])
    if (q.diff < 1) return(0.1)
    if (q.diff < 1.5) return(0.2)
    return(0.3)
}    

.getWellCoveredExons <- function(normal, tumor, min.coverage) {
    well.covered.exon.idx <- ((normal$average.coverage > min.coverage) & 
        (tumor$average.coverage > min.coverage)) | 
        ((normal$average.coverage > 1.5 * min.coverage) &  
        (tumor$average.coverage > 0.5 * min.coverage))
    #MR: fix for missing chrX/Y 
    well.covered.exon.idx[is.na(well.covered.exon.idx)] <- FALSE

    well.covered.exon.idx
}
        
# ExomeCNV version without the x11() calls 
.CNV.analyze2 <-
function(normal, tumor, logR=NULL, min.coverage=15, weights=NULL, sdundo=NULL,
undo.splits="sdundo", smooth=TRUE, alpha=0.01, sampleid=NULL, plot.cnv=TRUE,
max.segments=NULL, chr.hash=chr.hash, verbose=TRUE) {
    `%+%` <- function(x,y) paste(x,y,sep="")

    # first, do it for exons with enough coverage. MR: added less stringent 
    # cutoff in case normal looks great. these could be homozygous deletions 
    # in high purity samples
    well.covered.exon.idx <- .getWellCoveredExons(normal, tumor, 
        min.coverage)

    if (verbose) message("Removing ", sum(!well.covered.exon.idx), 
        " low coverage exons.")
    if (is.null(logR)) norm.log.ratio = calculateLogRatio(normal, tumor, verbose)
    else norm.log.ratio = logR
    
    if (is.null(sdundo)) {
        sdundo <- .getSDundo(norm.log.ratio[well.covered.exon.idx])
    }   
     
    CNA.obj <- CNA(norm.log.ratio[well.covered.exon.idx], 
        .strip.chr.name(normal$chr[well.covered.exon.idx], chr.hash), 
        (normal$probe_start[well.covered.exon.idx] + 
        normal$probe_end[well.covered.exon.idx])/2, data.type="logratio", 
        sampleid=sampleid)

    smoothed.CNA.obj = if (smooth) smooth.CNA(CNA.obj) else CNA.obj

    try.again <- 0

    while (try.again < 2) {
        if (verbose) message("Setting undo.SD parameter to ", sdundo)
        if (!is.null(weights)) { 
            weights <- weights[well.covered.exon.idx]
            # MR: this shouldn't happen. In doubt, count them as maximum 
            # (assuming that poorly performing exons are down-weighted)
            weights[is.na(weights)] <- max(weights, na.rm=TRUE)
            segment.smoothed.CNA.obj <- segment(smoothed.CNA.obj, 
                undo.splits=undo.splits, undo.SD=sdundo, 
                verbose=ifelse(verbose, 1, 0), alpha=alpha,weights=weights)
        } else {
            segment.smoothed.CNA.obj <- segment(smoothed.CNA.obj, 
                undo.splits=undo.splits, undo.SD=sdundo, 
                verbose=ifelse(verbose, 1, 0), alpha=alpha)
        } 
        if (is.null(max.segments) || nrow(segment.smoothed.CNA.obj$output) 
            < max.segments) break
        sdundo <- sdundo * 1.5
        try.again <- try.again + 1
    }

    if (plot.cnv) {
        plot(segment.smoothed.CNA.obj, plot.type="s")
        plot(segment.smoothed.CNA.obj, plot.type="w")
    }

    return(list(cna=segment.smoothed.CNA.obj, logR=norm.log.ratio))
}

.getSegSizes <- function(seg) {
    round(seg$loc.end-seg$loc.start+1)
}

