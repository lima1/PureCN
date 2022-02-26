#' Pool coverage from multiple samples
#'
#' Averages the coverage of a list of samples.
#'
#'
#' @param all.data List of normals, read with \code{\link{readCoverageFile}}.
#' @param remove.chrs Remove these chromosomes from the pool.
#' @param w \code{numeric(length(all.data))} vector of weights. If \code{NULL},
#' weight all samples equally.
#' @return A \code{data.frame} with the averaged coverage over all normals.
#' @author Markus Riester
#' @seealso \code{\link{readCoverageFile}}
#' @examples
#'
#' normal.coverage.file <- system.file("extdata", "example_normal.txt.gz",
#'     package = "PureCN")
#' normal2.coverage.file <- system.file("extdata", "example_normal2.txt.gz",
#'     package = "PureCN")
#' normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)
#' pool <- poolCoverage(lapply(normal.coverage.files, readCoverageFile),
#'      remove.chrs = c("chrX", "chrY"))
#'
#' @export poolCoverage
poolCoverage <- function(all.data, remove.chrs=c(), w = NULL) {
    pool <- all.data[[1]]

    if (length(all.data) == 1) {
        return(.removeChr(pool, remove.chrs))
    }
    if (is.null(w)) {
        w <- rep(1, length(all.data))
    } else if (length(w) != length(all.data)) {
        .stopUserError("all.data and w have different lengths.")
    }

    pool$coverage <- 0
    pool$counts <- 0

    for (i in seq_along(all.data)) {
        pool$coverage <- pool$coverage + (w[i] * all.data[[i]]$coverage)
        pool$counts <- pool$counts + (w[i] * all.data[[i]]$counts)
    }
    pool <- .addAverageCoverage(pool)
    return(.removeChr(pool, remove.chrs))
}

.removeChr <- function(pool, remove.chrs = c()) {
    idx <- seqnames(pool) %in% remove.chrs
    if (sum(idx)) {
        pool[idx]$coverage <- NA
        pool[idx]$average.coverage <- NA
    }
    pool
}
