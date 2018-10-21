#' PSCBS segmentation
#' 
#' Alternative segmentation function using the \code{PSCBS} package.  This
#' function is called via the \code{fun.segmentation} argument of
#' \code{\link{runAbsoluteCN}}.  The arguments are passed via
#' \code{args.segmentation}.
#' 
#' 
#' @param normal Coverage data for normal sample.
#' @param tumor Coverage data for tumor sample.
#' @param log.ratio Copy number log-ratios, one for each exon in coverage file.
#' @param seg If segmentation was provided by the user, this data structure
#' will contain this segmentation. Useful for minimal segmentation functions.
#' Otherwise PureCN will re-segment the data. This segmentation function
#' ignores this user provided segmentation.
#' @param plot.cnv Segmentation plots.
#' @param sampleid Sample id, used in output files.
#' @param interval.weight.file Can be used to assign weights to intervals. NOT
#' SUPPORTED YET in segmentation. Will remove targets with weight below 1/3.
#' @param target.weight.file Deprecated.
#' @param alpha Alpha value for CBS, see documentation for the \code{segment}
#' function.
#' @param undo.SD \code{undo.SD} for CBS, see documentation of the
#' \code{segment} function. If \code{NULL}, try to find a sensible default.
#' @param flavor Flavor value for PSCBS. See \code{segmentByNonPairedPSCBS}.
#' @param tauA tauA argument for PSCBS. See \code{segmentByNonPairedPSCBS}.
#' @param vcf Optional VCF object with germline allelic ratios.
#' @param tumor.id.in.vcf Id of tumor in case multiple samples are stored in
#' VCF.
#' @param normal.id.in.vcf Id of normal in in VCF. If \code{NULL}, 
#' use unpaired PSCBS.
#' @param max.segments If not \code{NULL}, try a higher \code{undo.SD} 
#' parameter if number of segments exceeds the threshold.
#' @param prune.hclust.h Height in the \code{hclust} pruning step. Increasing
#' this value will merge segments more aggressively. If \code{NULL}, try to 
#' find a sensible default.
#' @param prune.hclust.method Cluster method used in the \code{hclust} pruning
#' step. See documentation for the \code{hclust} function.
#' @param chr.hash Mapping of non-numerical chromsome names to numerical names
#' (e.g. chr1 to 1, chr2 to 2, etc.). If \code{NULL}, assume chromsomes are 
#' properly ordered.
#' @param centromeres A \code{GRanges} with centromere positions.
#' If not \code{NULL}, add breakpoints at centromeres. 
#' @param \dots Additional parameters passed to the 
#' \code{segmentByNonPairedPSCBS} function.
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
#' Olshen et al. (2011). Parent-specific copy number in paired tumor-normal 
#' studies using circular binary segmentation. Bioinformatics.
#' @seealso \code{\link{runAbsoluteCN}}
#' @examples
#' 
#' normal.coverage.file <- system.file("extdata", "example_normal_tiny.txt", 
#'     package="PureCN")
#' tumor.coverage.file <- system.file("extdata", "example_tumor_tiny.txt", 
#'     package="PureCN")
#' vcf.file <- system.file("extdata", "example.vcf.gz",
#'     package="PureCN")
#' 
#' # The max.candidate.solutions, max.ploidy and test.purity parameters are set to
#' # non-default values to speed-up this example.  This is not a good idea for real
#' # samples.
#'  ret <-runAbsoluteCN(normal.coverage.file=normal.coverage.file, 
#'      tumor.coverage.file=tumor.coverage.file, vcf.file=vcf.file, 
#'      sampleid="Sample1",  genome="hg19",
#'      fun.segmentation=segmentationPSCBS, max.ploidy=4,
#'      test.purity=seq(0.3,0.7,by=0.05), max.candidate.solutions=1)
#' 
#' @export segmentationPSCBS
segmentationPSCBS <- function(normal, tumor, log.ratio, seg, plot.cnv, 
    sampleid, interval.weight.file = NULL, target.weight.file = NULL, alpha = 0.005, 
    undo.SD = NULL, flavor = "tcn&dh", tauA = 0.03, vcf = NULL,
    tumor.id.in.vcf = 1, normal.id.in.vcf = NULL, max.segments = NULL,
    prune.hclust.h = NULL, prune.hclust.method = "ward.D", chr.hash = NULL,
    centromeres = NULL, ...) {

    if (!requireNamespace("PSCBS", quietly = TRUE)) {
        .stopUserError("segmentationPSCBS requires the PSCBS package.")
    }

    # TODO remove in 1.14.0
    if (is.null(interval.weight.file) && !is.null(target.weight.file)) {
        interval.weight.file <- target.weight.file
        flog.warn("target.weight.file was renamed to interval.weight.file.")
    }

    if (is.null(chr.hash)) chr.hash <- .getChrHash(seqlevels(tumor))

    interval.weights <- NULL
    well.covered.exon.idx <- rep(TRUE, length(tumor))
    if (!is.null(interval.weight.file)) {
        interval.weights <- read.delim(interval.weight.file, as.is=TRUE)
        interval.weights <- interval.weights[match(as.character(tumor), 
            interval.weights[,1]),2]
         flog.info("Interval weights found, but currently not supported by PSCBS. %s",
            "Will simply exclude intervals with low weight.")
        lowWeightIntervals <- interval.weights < 1/3
        well.covered.exon.idx[which(lowWeightIntervals)] <- FALSE
    }

    #MR: fix for missing chrX/Y 
    well.covered.exon.idx[is.na(well.covered.exon.idx)] <- FALSE
    tumor <- tumor[well.covered.exon.idx]
    log.ratio <- log.ratio[well.covered.exon.idx]
    ov <- findOverlaps(vcf, tumor)
    d.f <- cbind(as.data.frame(tumor[subjectHits(ov)]), 
        CT=2 ^ (log.ratio+1)[subjectHits(ov)], 
        betaT=unlist(geno(vcf[queryHits(ov)])$FA[,tumor.id.in.vcf]), 
        betaN=NA,
        x=start(vcf[queryHits(ov)]))
    
    if (!is.null(normal.id.in.vcf)) {
        d.f$betaN <- unlist(geno(vcf[queryHits(ov)])$FA[,normal.id.in.vcf])
    }
         
    d.f.2 <- cbind(as.data.frame(tumor[-subjectHits(ov)]), 
        CT=2 ^ (log.ratio+1)[-subjectHits(ov)], betaT=NA, betaN=NA,
        x=start(tumor[-subjectHits(ov)]))
    
    d.f <- rbind(d.f, d.f.2)
    colnames(d.f)[1] <- "chromosome"
    d.f <- d.f[order(.strip.chr.name(d.f[,1], chr.hash), d.f$x),]
    d.f$chromosome <- .strip.chr.name(d.f$chromosome, chr.hash)

    if (is.null(undo.SD)) {
        undo.SD <- .getSDundo(log.ratio)
        flog.info("Setting undo.SD parameter to %f.", undo.SD)
    }   
    knownSegments <- .PSCBSgetKnownSegments(centromeres, chr.hash)
    seg <- PSCBS::segmentByNonPairedPSCBS(d.f, tauA=tauA, 
        flavor=flavor, undoTCN=undo.SD, knownSegments=knownSegments, 
        min.width=3,alphaTCN=alpha, ...)

    if (plot.cnv) PSCBS::plotTracks(seg)
    seg <- .PSCBSoutput2DNAcopy(seg, sampleid)

    if (!is.null(vcf)) {
        seg <- .pruneByHclust(seg, vcf, tumor.id.in.vcf, h=prune.hclust.h, 
            method=prune.hclust.method, chr.hash=chr.hash)
    }
    seg
}

.PSCBSgetKnownSegments <- function(centromeres, chr.hash) {
    if (is.null(centromeres)) return(NULL)
    knownSegments <- data.frame(centromeres)
    colnames(knownSegments)[1] <- "chromosome"
    knownSegments$length <- knownSegments$end-knownSegments$start+1
    knownSegments$chromosome <- .strip.chr.name(knownSegments$chromosome,
        chr.hash)
    PSCBS::gapsToSegments(knownSegments)
}

.PSCBSoutput2DNAcopy <- function(seg, sampleid) {
    sx <- cbind(ID=sampleid, seg$output[!is.na(seg$output$tcnMean),])
    sx <- sx[,c("ID", "chromosome", "tcnStart", "tcnEnd", "tcnNbrOfLoci", 
        "tcnMean")]
    colnames(sx) <- c("ID", "chrom", "loc.start",  "loc.end", "num.mark", 
        "seg.mean")
    sx$seg.mean <- log2(sx$seg.mean/2)
    sx
}
