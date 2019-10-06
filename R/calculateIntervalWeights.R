#' Calculate interval weights
#' 
#' Creates an interval weight file useful for segmentation. Requires a set of 
#' coverage files from normal samples. Interval weights will be
#' set proportional to the inverse of coverage standard deviation across all
#' normals. Intervals with high variance in coverage in the pool of normals are
#' thus down-weighted.
#' 
#' This function is now automatically called by \code{\link{createNormalDatabase}}
#' and is thus deprecated.
#' 
#' @param normalDB Database of normal samples, created with
#' \code{\link{createNormalDatabase}}.
#' @param interval.weight.file Output filename.
#' @param top.quantile Cap weight at the specified quantile. Intervals
#' with standard deviation smaller than this value won't have a higher weight
#' than intervals at this quantile.
#' @param plot Diagnostics plot, useful to tune parameters.
#' @param normal.coverage.files Deprecated.
#' @return A normalDB object with following slots added
#' \item{sd$log.ratios}{\code{GRanges} with all log.ratios.}
#' \item{sd$weights}{\code{GRanges} with interval weights.}
#' @author Markus Riester
#' 
#' @export calculateIntervalWeights
calculateIntervalWeights <- function(normalDB,
interval.weight.file = NULL, top.quantile = 0.7, plot = FALSE, 
normal.coverage.files = NULL) {
    .Deprecated("createNormalDatabase")

    # TODO, defunct in 1.18
    old_method <- FALSE
    if (!is.null(normal.coverage.files) && missing(normalDB)) {
        flog.warn("normal.coverage.files is deprecated, provide normalDB instead.")
        old_method = TRUE
    } else if (class(normalDB) == "character" && all(sapply(normalDB, file.exists)))  {
        flog.warn("Providing normal coverage files is deprecated. Provide a normalDB instead.")
        normal.coverage.files <- normalDB
        old_method = TRUE
    } else {
        normal.coverage.files <- normalDB[["normal.coverage.files"]]
    }    
    flog.info("Loading coverage data...")
    normal.coverage <- lapply(normal.coverage.files,  readCoverageFile)
    .calculateIntervalWeights(normalDB, normal.coverage, interval.weight.file, top.quantile,
        plot, old_method) 
}    

.calculateIntervalWeights <- function(normalDB, normal.coverage, 
interval.weight.file = NULL, top.quantile = 0.7, plot = FALSE, 
old_method = FALSE) {
    
    tumor.coverage <- list(poolCoverage(normal.coverage, 
        w = rep(1, length(normal.coverage)) / length(normal.coverage)))

    if (old_method) {
        lrs <- lapply(tumor.coverage, function(tc) sapply(normal.coverage, 
                function(nc) calculateLogRatio(nc, tc)))
    } else {
        lrs <- lapply(normal.coverage, function(x) calculateTangentNormal(x, normalDB)$log.ratio)
    }

    lrs <- do.call(cbind, lrs)

    lrs[is.infinite(lrs)] <- NA

    intervals <- normal.coverage[[1]]
    mcols(intervals) <- NULL

    lrs.sd <- apply(lrs, 1, sd, na.rm = TRUE)
    lrs.cnt.na <- apply(lrs, 1, function(x) sum(is.na(x)))
    # get the top.quantile % of sd by chromosome and use this to normalize weight=1
    chrom <-  as.character(seqnames(intervals))
    sdCutoffByChr <- sapply(split(lrs.sd, chrom), quantile, probs = top.quantile, 
        names = FALSE, na.rm = TRUE)[chrom]

    zz <- sdCutoffByChr / lrs.sd
    zz[zz > 1] <- 1
    idx <- is.na(zz) | lrs.cnt.na > ncol(lrs) / 3
    zz[idx] <- min(zz, na.rm = TRUE)

    ret <- list(
                log.ratios = GRanges(intervals,,, DataFrame(lrs)),
                weights = GRanges(intervals,,, DataFrame(weights = zz)))
    
    if (!is.null(interval.weight.file)) {
        ret_output <- data.frame(
            Target = as.character(ret$weights), 
            weights = ret$weights$weights)

        fwrite(ret_output, file = interval.weight.file, row.names = FALSE, 
                    quote = FALSE, sep = "\t")
    }
    if (plot) .plotIntervalWeights(lrs.sd, width(tumor.coverage[[1]]), 
        tumor.coverage[[1]]$on.target)
    if (old_method) return(NULL)    
    normalDB$sd <- ret
    normalDB
}

#' Calculate target weights
#' 
#' This function is defunct, use \code{\link{calculateIntervalWeights}}
#' instead.
#'
#' @author Markus Riester
#'
#' @export createTargetWeights
createTargetWeights <- function() {
    .Defunct("calculateIntervalWeights")
}

.plotIntervalWeights <- function(lrs.sd, width, on.target) {
    par(mfrow = c(1, 2))
    plot(width[on.target], lrs.sd[on.target], ylim = c(0,2),
        xlab = "Interval Width", ylab = "log2 ratio sd.", main = "On-Target")
    if (sum(!on.target)) {
        plot(width[!on.target], lrs.sd[!on.target], col = "red", 
             ylim = c(0, 2), xlab = "Interval Width", ylab = "log2 ratio sd.",
             main = "Off-Target")
    }
}
