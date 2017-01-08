#' PSCBS segmentation
#' 
#' Alternative segmentation function using the \code{PSCBS} package.  This
#' function is called via the \code{fun.segmentation} argument of
#' \code{\link{runAbsoluteCN}}.  The arguments are passed via
#' \code{args.segmentation}.
#' 
#' 
#' @param normal GATK coverage data for normal sample.
#' @param tumor GATK coverage data for tumor sample.
#' @param log.ratio Copy number log-ratios, one for each exon in coverage file.
#' @param seg If segmentation was provided by the user, this data structure
#' will contain this segmentation. Useful for minimal segmentation functions.
#' Otherwise PureCN will re-segment the data. This segmentation function
#' ignores this user provided segmentation.
#' @param plot.cnv Segmentation plots.
#' @param min.coverage Minimum coverage in both normal and tumor.
#' @param sampleid Sample id, used in output files.
#' @param target.weight.file Can be used to assign weights to targets. NOT
#' SUPPORTED YET in segmentation. Will remove targets with weight below 1/3.
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
#' @param centromeres A \code{data.frame} with centromere positions in first
#' three columns.  If not \code{NULL}, add breakpoints at centromeres. 
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
#' normal.coverage.file <- system.file("extdata", "example_normal.txt", 
#'     package="PureCN")
#' tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
#'     package="PureCN")
#' vcf.file <- system.file("extdata", "example_vcf.vcf", 
#'     package="PureCN")
#' gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
#'     package="PureCN")
#' 
#' # The max.candidate.solutions, max.ploidy and test.purity parameters are set to
#' # non-default values to speed-up this example.  This is not a good idea for real
#' # samples.
#'  ret <-runAbsoluteCN(normal.coverage.file=normal.coverage.file, 
#'      tumor.coverage.file=tumor.coverage.file, vcf.file=vcf.file, genome="hg19",
#'      sampleid="Sample1", gc.gene.file=gc.gene.file,
#'      fun.segmentation=segmentationPSCBS, max.ploidy=4,
#'      test.purity=seq(0.3,0.7,by=0.05), max.candidate.solutions=1)
#' 
#' @export segmentationPSCBS
segmentationPSCBS <- function(normal, tumor, log.ratio, seg, plot.cnv, 
    min.coverage, sampleid, target.weight.file = NULL, alpha = 0.005, undo.SD =
    NULL, flavor = "tcn&dh", tauA = 0.03, vcf = NULL,
    tumor.id.in.vcf = 1, normal.id.in.vcf = NULL, max.segments = NULL,
    prune.hclust.h = NULL, prune.hclust.method = "ward.D", chr.hash = NULL,
    centromeres = NULL, ...) {

    debug <- TRUE
        
    if (!requireNamespace("PSCBS", quietly = TRUE)) {
        .stopUserError("segmentationPSCBS requires the PSCBS package.")
    }

    if (is.null(chr.hash)) chr.hash <- .getChrHash(tumor$chr)

    well.covered.exon.idx <- ((normal$average.coverage > min.coverage) &
        (tumor$average.coverage > min.coverage)) | 
        ((normal$average.coverage > 1.5 * min.coverage) &  
        (tumor$average.coverage > 0.5 * min.coverage))

    target.weights <- NULL
    if (!is.null(target.weight.file)) {
        target.weights <- read.delim(target.weight.file, as.is=TRUE)
        target.weights <- target.weights[match(as.character(tumor[,1]), 
            target.weights[,1]),2]
         flog.info("Target weights found, but currently not supported by PSCBS. %s",
            "Will simply exclude targets with low weight.")
        lowWeightTargets <- target.weights < 1/3
        well.covered.exon.idx[which(lowWeightTargets)] <- FALSE
    }

    #MR: fix for missing chrX/Y 
    well.covered.exon.idx[is.na(well.covered.exon.idx)] <- FALSE
    tumor <- tumor[well.covered.exon.idx,]
    log.ratio <- log.ratio[well.covered.exon.idx]
    exon.gr <- GRanges(seqnames=tumor$chr, 
        IRanges(start=tumor$probe_start, end=tumor$probe_end))
    ov <- findOverlaps(vcf, exon.gr)
    d.f <- cbind(tumor[subjectHits(ov),], 
        CT=2^(log.ratio+1)[subjectHits(ov)], 
        betaT=unlist(geno(vcf[queryHits(ov)])$FA[,tumor.id.in.vcf]), 
        betaN=NA,
        x=start(vcf[queryHits(ov)]) )
    
    if (!is.null(normal.id.in.vcf)) {
        d.f$betaN <- unlist(geno(vcf[queryHits(ov)])$FA[,normal.id.in.vcf])
    }
         
    d.f.2 <- cbind(tumor[-subjectHits(ov),], 
        CT=2^(log.ratio+1)[-subjectHits(ov)], betaT=NA, betaN=NA,
        x=tumor$probe_start[-subjectHits(ov)] )
    
    d.f.3 <- rbind(d.f, d.f.2)
    d.f.3 <- d.f.3[order(.strip.chr.name(d.f.3$chr, chr.hash), d.f.3$x),]
    d.f <- d.f.3
    colnames(d.f)[2] <- "chromosome"
    d.f$chromosome <- .strip.chr.name(d.f$chromosome, chr.hash)
    #if (!is.null(normal.id.in.vcf)) {
    #    seg <- PSCBS::segmentByPairedPSCBS(d.f, tauA=tauA, 
    #        flavor=flavor, ...)
    #} else {
    if (is.null(undo.SD)) {
        undo.SD <- .getSDundo(log.ratio)
        flog.info("Setting undo.SD parameter to %f.", undo.SD)
    }   
    knownSegments <- NULL
    if (!is.null(centromeres)) {
        knownSegments <- centromeres
        colnames(knownSegments) <- c("chromosome", "start", "end")
        knownSegments$length <- knownSegments$end-knownSegments$start+1
        knownSegments$chromosome <- .strip.chr.name(knownSegments$chromosome,
            chr.hash)
        knownSegments <- PSCBS::gapsToSegments(knownSegments)
    }    
    seg <- PSCBS::segmentByNonPairedPSCBS(d.f, tauA=tauA, 
        flavor=flavor, undoTCN=undo.SD, knownSegments=knownSegments, 
        min.width=3,alphaTCN=alpha, ...)
    #}    
    if (plot.cnv) PSCBS::plotTracks(seg)
    x <- .PSCBSoutput2DNAcopy(seg, sampleid)

    if (!is.null(vcf)) {
        x <- .pruneByHclust(x, vcf, tumor.id.in.vcf, h=prune.hclust.h, 
            method=prune.hclust.method, chr.hash=chr.hash)
    }
    x$cna$output
}

.PSCBSoutput2DNAcopy <- function(seg, sampleid) {
    sx <- cbind(ID=sampleid, seg$output[!is.na(seg$output$tcnMean),])
    sx <- sx[,c("ID", "chromosome", "tcnStart", "tcnEnd", "tcnNbrOfLoci", 
        "tcnMean")]
    colnames(sx) <- c("ID", "chrom", "loc.start",  "loc.end", "num.mark", 
        "seg.mean")
    sx$seg.mean <- log2(sx$seg.mean/2)
    list(cna=list(output=sx))
}
