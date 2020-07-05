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
#' @param genome Genome version, for example hg19. Needed to get centromere
#' positions.
#' @param plot.cnv Segmentation plots.
#' @param min.interval.weight Can be used to ignore intervals with low weights.
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
#' @references Nilsen G., Liestol K., Van Loo P., Vollan H., Eide M., Rueda O.,
#' Chin S., Russell R., Baumbusch L., Caldas C., Borresen-Dale A., 
#' Lingjaerde O. (2012). "Copynumber: Efficient algorithms for single- and
#' multi-track copy number segmentation." BMC Genomics, 13(1), 591.
#'
#' @seealso \code{\link{runAbsoluteCN}}
#' @examples
#' 
#' normal1.coverage.file <- system.file("extdata", "example_normal.txt", 
#'     package="PureCN")
#' normal2.coverage.file <- system.file("extdata", "example_normal2.txt", 
#'     package="PureCN")
#' tumor1.coverage.file <- system.file("extdata", "example_tumor.txt",
#'     package="PureCN")
#' tumor2.coverage.file <- system.file("extdata", "example_tumor2.txt",
#'     package="PureCN")
#'
#' normal.coverage.files <- c(normal1.coverage.file, normal2.coverage.file)
#' tumor.coverage.files <- c(tumor1.coverage.file, tumor2.coverage.file)
#'
#' normalDB <- createNormalDatabase(normal.coverage.files)
#'
#' seg <- processMultipleSamples(tumor.coverage.files, 
#'          sampleids = c("Sample1", "Sample2"),
#'          normalDB = normalDB,
#'          genome = "hg19")
#'
#' @export processMultipleSamples
processMultipleSamples <- function(tumor.coverage.files, sampleids, normalDB, 
    num.eigen = 20, genome, plot.cnv = TRUE, w = NULL,
    min.interval.weight = 1/3,
    max.segments = NULL, chr.hash = NULL, centromeres = NULL, ...) {

    if (!requireNamespace("copynumber", quietly = TRUE)) {
        .stopUserError("processMultipleSamples requires the copynumber package.")
    }
    tumor.coverage.files <- normalizePath(tumor.coverage.files)
    tumors <- lapply(.readNormals(tumor.coverage.files), 
        calculateTangentNormal, normalDB, num.eigen = num.eigen)

    interval.weights <- NULL
    intervalsUsed <- rep(TRUE, length(tumors[[1]]))
    if (!is.null(normalDB$sd$weights) && !is.null(min.interval.weight)) {
       interval.weights <- normalDB$sd$weights$weights
       if (length(interval.weights) != length(tumors[[1]])) {
            # should not happen because it is checked upstream
            .stopRuntimeError("interval.weights and tumors does not align.")
       }    
       flog.info("Interval weights found, but currently not supported by copynumber. %s",
            "Will simply exclude intervals with low weight.")
       lowWeightIntervals <- interval.weights < min.interval.weight
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
    # ignore arms with only 2 or fewer probes
    carms <- paste0(lrs[,1], arms)
    intervalsUsed <- as.logical(intervalsUsed & complete.cases(lrs) 
                                & !is.na(arms)) & table(carms)[carms] > 2

    lrs <- lrs[intervalsUsed,]
    arms <- arms[intervalsUsed]
    lrsw <- copynumber::winsorize(lrs, arms = arms, verbose = FALSE)
    if (is.null(w)) {
        w <- 1
        dupr <- vapply(tumors, function(x) 
                       median(x[x$on.target]$duplication.rate, na.rm = TRUE),
                       double(1))
        if (!sum(is.na(dupr)) && min(dupr, na.rm = TRUE) > 0) { 
            w <- (1 / dupr)
            w <- w / max(w)

            flog.info("Setting weights by duplication rate. Lowest weight for %s (%.2f), heighest for %s.",
                sampleids[which.min(w)], min(w), sampleids[which.max(w)])
        }
    }
    lrsm <- copynumber::multipcf(lrsw, arms = arms, w = w, ...)
    if (plot.cnv) {
        copynumber::plotHeatmap(segments = lrsm, upper.lim = 1)
        copynumber::plotGenome(lrsw, segments = lrsm, onefile = TRUE)
    }
    idx.enough.markers <- lrsm$n.probes > 1
    rownames(lrsm) <- NULL
    lrsm[idx.enough.markers,]
    #transform to DNAcopy format
    m <- data.table::melt(data.table(lrsm), id.vars=1:5)
    m <- m[, c(6,1,3,4,5,7)]
    colnames(m) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
    data.frame(m)
}
