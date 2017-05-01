#' Calculate coverage log-ratio of tumor vs. normal
#' 
#' This function is automatically called by \code{\link{runAbsoluteCN}} when
#' normal and tumor coverage are provided (and not a segmentation file or
#' target-level log-ratios). This function is therefore normally not called by
#' the user.
#' 
#' 
#' @param normal Normal coverage read in by the \code{\link{readCoverageFile}}
#' function.
#' @param tumor Tumor coverage read in by the \code{\link{readCoverageFile}}
#' function.
#' @return \code{numeric(length(tumor))}, tumor vs. normal copy number log-ratios
#' for all targets.
#' @author Markus Riester
#' @examples
#' 
#' normal.coverage.file <- system.file("extdata", "example_normal.txt", 
#'     package="PureCN")
#' tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
#'     package="PureCN")
#' normal <- readCoverageFile(normal.coverage.file)
#' tumor <- readCoverageFile(tumor.coverage.file)
#' log.ratio <- calculateLogRatio(normal, tumor)
#' 
#' @export calculateLogRatio
calculateLogRatio <- function(normal, tumor) {
    # make sure that normal and tumor align
    if (!identical(as.character(normal), as.character(tumor))) {
        .stopUserError("Interval files in normal and tumor different.")
    }
    flog.info("Mean coverages: %.0fX (tumor) %.0fX (normal).", 
        mean(tumor$average.coverage, na.rm=TRUE), 
        mean(normal$average.coverage, na.rm=TRUE))
    total.cov.normal <- sum(as.numeric(normal$coverage), na.rm = TRUE)
    total.cov.tumor <- sum(as.numeric(tumor$coverage), na.rm = TRUE)

    log.ratio <- log2(tumor$average.coverage/normal$average.coverage) + 
                 log2(total.cov.normal/total.cov.tumor)

    mean.log.ratio <- mean(subset(log.ratio, !is.infinite(log.ratio)), 
        na.rm = TRUE)
    # calibrate
    log.ratio <- log.ratio - mean.log.ratio
    log.ratio
}
