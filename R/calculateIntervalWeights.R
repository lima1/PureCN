#' Calculate interval weights
#' 
#' Creates an interval weight file useful for segmentation. Requires a set of 
#' coverage files from normal samples. Interval weights will be
#' set proportional to the inverse of coverage standard deviation across all
#' normals. Intervals with high variance in coverage in the pool of normals are
#' thus down-weighted.
#' 
#' 
#' @param normal.coverage.files A set of normal coverage samples
#' to estimate target log-ratio standard deviations. 
#' @param interval.weight.file Output filename.
#' @param plot Diagnostics plot, useful to tune parameters.
#' @return A \code{data.frame} with interval weights.
#' @author Markus Riester
#' @examples
#' 
#' interval.weight.file <- "interval_weights.txt"
#' normal.coverage.file <- system.file("extdata", "example_normal.txt", 
#'     package="PureCN")
#' normal2.coverage.file <- system.file("extdata", "example_normal2.txt", 
#'     package="PureCN")
#' normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)
#' 
#' calculateIntervalWeights(normal.coverage.files, interval.weight.file)
#' 
#' @export calculateIntervalWeights
calculateIntervalWeights <- function(normal.coverage.files,
interval.weight.file, plot = FALSE) {
    flog.info("Loading coverage data...")
    normal.coverage <- lapply(normal.coverage.files,  readCoverageFile)
    
    tumor.coverage <- list(poolCoverage(normal.coverage, 
        w = rep(1, length(normal.coverage)) / length(normal.coverage)))

    lrs <- lapply(tumor.coverage, function(tc) sapply(normal.coverage, 
            function(nc) calculateLogRatio(nc, tc)))

    lrs <- do.call(cbind, lrs)

    lrs[is.infinite(lrs)] <- NA

    lrs.sd <- apply(lrs, 1, sd, na.rm = TRUE)
    lrs.cnt.na <- apply(lrs, 1, function(x) sum(is.na(x)))
    # get the 70% of sd by chromosome and use this to normalize weight=1
    chrom <-  as.character(seqnames(tumor.coverage[[1]]))
    sdCutoffByChr <- sapply(split(lrs.sd, chrom), quantile, probs = 0.7, 
        names = FALSE, na.rm = TRUE)[chrom]

    zz <- sdCutoffByChr / lrs.sd
    zz[zz > 1] <- 1
    idx <- is.na(zz) | lrs.cnt.na > ncol(lrs) / 3
    zz[idx] <- min(zz, na.rm = TRUE)
    ret <- data.frame(Target = as.character(tumor.coverage[[1]]), Weights = zz)
    
    write.table(ret, file = interval.weight.file, row.names = FALSE, 
                quote = FALSE, sep = "\t")
    if (plot) .plotIntervalWeights(lrs.sd, width(tumor.coverage[[1]]), 
        tumor.coverage[[1]]$on.target)
    invisible(ret)
}

#' Calculate target weights
#' 
#' This function is deprecated, use \code{\link{calculateIntervalWeights}}
#' instead.
#'
#' @param normal.coverage.files A set of normal coverage samples
#' to estimate target log-ratio standard deviations. 
#' @param target.weight.file Output filename.
#' @param plot Diagnostics plot, useful to tune parameters.
#' @return A \code{data.frame} with target weights.
#' @author Markus Riester
#'
#' @export createTargetWeights
createTargetWeights <- function(normal.coverage.files,
target.weight.file, plot = FALSE) {
    .Deprecated("calculateIntervalWeights")
    calculateIntervalWeights(normal.coverage.files, target.weight.file, plot)
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
