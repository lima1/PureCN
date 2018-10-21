#' CBS segmentation
#' 
#' The default segmentation function. This function is called via the
#' \code{fun.segmentation} argument of \code{\link{runAbsoluteCN}}.  The
#' arguments are passed via \code{args.segmentation}.
#' 
#' 
#' @param seg If segmentation was provided by the user, this data structure
#' will contain this segmentation. Useful for minimal segmentation functions.
#' Otherwise PureCN will re-segment the data. This segmentation function
#' ignores this user provided segmentation.
#' @param vcf Optional \code{CollapsedVCF} object with germline allelic ratios.
#' @param tumor.id.in.vcf Id of tumor in case multiple samples are stored in
#' VCF.
#' @param normal.id.in.vcf Id of normal in in VCF. Currently not used.
#' @param prune.hclust.h Height in the \code{hclust} pruning step. Increasing
#' this value will merge segments more aggressively. If NULL, try to find a
#' sensible default.
#' @param prune.hclust.method Cluster method used in the \code{hclust} pruning
#' step. See documentation for the \code{hclust} function.
#' @param chr.hash Mapping of non-numerical chromsome names to numerical names
#' (e.g. chr1 to 1, chr2 to 2, etc.). If \code{NULL}, assume chromsomes are
#' properly ordered.
#' @param ... Currently unused arguments provided to other segmentation
#' functions.
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
#' normal.coverage.file <- system.file("extdata", "example_normal_tiny.txt", 
#'     package="PureCN")
#' tumor.coverage.file <- system.file("extdata", "example_tumor_tiny.txt", 
#'     package="PureCN")
#' vcf.file <- system.file("extdata", "example.vcf.gz", 
#'     package="PureCN")
#' interval.file <- system.file("extdata", "example_intervals_tiny.txt", 
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
#' @export segmentationHclust
segmentationHclust <- function(seg, 
    vcf = NULL, tumor.id.in.vcf = 1, normal.id.in.vcf = NULL,
    prune.hclust.h = NULL, prune.hclust.method = "ward.D",
    chr.hash = NULL, ...) {
    
    if (!is.null(vcf)) {
        if (is.null(chr.hash)) chr.hash <- .getChrHash(seqlevels(vcf))
        seg <- .pruneByHclust(seg, vcf, tumor.id.in.vcf, 
            h = prune.hclust.h, 
            method=prune.hclust.method, chr.hash = chr.hash)
    }
    idx.enough.markers <- seg$num.mark > 1
    rownames(seg) <- NULL
    seg[idx.enough.markers,]
}
