#' Multi sample normalization and segmentation
#' 
#' This function performs normalization and segmentation when multiple 
#' for the same patient are available. 
#' 
#' 
#' @param tumor.coverage.files Coverage data for tumor samples.
#' @param sampleids Sample ids, used in output files.
#' @param normalDB Database of normal samples, created with
#' \code{\link{createNormalDatabase}}.
#' @param num.eigen Number of eigen vectors used.
#' @param genome Genome version, for example hg19. See \code{readVcf}.
#' @param plot.cnv Segmentation plots.
#' @param interval.weight.file Can be used to assign weights to intervals.
#' @param w Weight of samples. Can be used to downweight poor quality samples.
#' If \code{NULL}, sets to inverse of median on-target duplication rate if
#' available, otherwise does not do any weighting.
#' @param max.segments If not \code{NULL}, try a higher \code{undo.SD}
#' parameter if number of segments exceeds the threshold.
#' @param chr.hash Mapping of non-numerical chromsome names to numerical names
#' (e.g. chr1 to 1, chr2 to 2, etc.). If \code{NULL}, assume chromsomes are
#' properly ordered.
#' @param centromeres A \code{GRanges} object with centromere positions.
#' @param ... Arguments passed to the segmentation function.
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
#' @export processMultipleSamples
processMultipleSamples <- function(tumor.coverage.files, sampleids, normalDB, 
    num.eigen = 20, genome, plot.cnv = TRUE, w = NULL,
    interval.weight.file = NULL, max.segments = NULL, 
    chr.hash = NULL, centromeres = NULL, ...) {

    if (!requireNamespace("copynumber", quietly = TRUE)) {
        .stopUserError("processMultipleSamples requires the copynumber package.")
    }
    tumor.coverage.files <- normalizePath(tumor.coverage.files)
    tumors <- lapply(.readNormals(tumor.coverage.files), 
        calculateTangentNormal, normalDB, num.eigen = num.eigen)

    interval.weights <- NULL
    intervalsUsed <- rep(TRUE, length(tumors[[1]]))
    if (!is.null(interval.weight.file)) {
        interval.weights <- read.delim(interval.weight.file, as.is=TRUE)
        interval.weights <- interval.weights[match(as.character(tumors[[1]]), 
            interval.weights[,1]),2]
         flog.info("Interval weights found, but currently not supported by copynumber. %s",
            "Will simply exclude intervals with low weight.")
        lowWeightIntervals <- interval.weights < 1/3
        intervalsUsed[which(lowWeightIntervals)] <- FALSE
    }

    #MR: fix for missing chrX/Y 
    intervalsUsed[is.na(intervalsUsed)] <- FALSE
    if (is.null(chr.hash)) chr.hash <- .getChrHash(seqlevels(tumors[[1]]))
    intervalsUsed <- .filterIntervalsChrHash(intervalsUsed, tumors[[1]], chr.hash)
    centromeres <- .getCentromerePositions(centromeres, genome, 
            if (is.null(tumors[[1]])) NULL else seqlevelsStyle(tumors[[1]]))

    armLocations <- .getArmLocations(tumors[[1]], chr.hash, centromeres)
    armLocationsGr <- GRanges(armLocations)
    arms <- armLocationsGr$arm[findOverlaps(tumors[[1]], armLocationsGr, select="first")]
    lrs <- data.frame(do.call(cbind, lapply(tumors, function(x) unlist(x$log.ratio))))
    colnames(lrs) <- sampleids
    lrs <- data.frame(chrom = .strip.chr.name(seqnames(tumors[[1]]), chr.hash), 
                      pos = start(tumors[[1]]), lrs)
    intervalsUsed <- as.logical(intervalsUsed & complete.cases(lrs) & !is.na(arms))

    lrs <- lrs[intervalsUsed,]
    arms <- arms[intervalsUsed]
    lrsw <- copynumber::winsorize(lrs, arms = arms)
    if (is.null(w)) {
        w <- 1
        dupr <- sapply(tumors, function(x) median(x[x$on.target]$duplication.rate, na.rm = TRUE))
        if (!sum(is.na(dupr)) && min(dupr, na.rm = TRUE) > 0) { 
            w <- (1/dupr)
            w <- w/max(w)

            flog.info("Setting weights by duplication rate. Lowest weight for %s (%.2f), heighest for %s.",
                sampleids[which.min(w)], min(w), sampleids[which.max(w)])
        }
    }
    lrsm <- copynumber::multipcf(lrsw, arms = arms, w = w, ...)
    if (plot.cnv) {
        copynumber::plotGenome(lrsw, segments = lrsm, onefile = TRUE)
        copynumber::plotHeatmap(segments = lrsm, upper.lim = 1)
    }
    idx.enough.markers <- lrsm$n.probes > 1
    rownames(lrsm) <- NULL
    lrsm[idx.enough.markers,]
    #transform to DNAcopy format
    m <- data.table::melt(lrsm, id.vars=1:5)
    m <- m[, c(6,1,3,4,5,7)]
    colnames(m) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
    m
}

