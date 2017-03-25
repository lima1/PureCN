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
#' normal.coverage.file <- system.file("extdata", "example_normal.txt", 
#'     package="PureCN")
#' normal2.coverage.file <- system.file("extdata", "example_normal2.txt", 
#'     package="PureCN")
#' normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)
#' tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
#'     package="PureCN")
#' 
#' normalDB <- createNormalDatabase(normal.coverage.files)
#' 
#' # get the best 2 normals and average them
#' best.normal.coverage.files <- findBestNormal(tumor.coverage.file, normalDB, 
#'     num.normals=2)
#' pool <- poolCoverage(lapply(best.normal.coverage.files, readCoverageFile),
#'      remove.chrs=c("chrX", "chrY"))
#' 
#' @export poolCoverage
poolCoverage <- function(all.data, remove.chrs=c(), w = NULL) { 
    pool = all.data[[1]]
    if (length(all.data) == 1) {
        return(.removeChr(pool, remove.chrs))
    }
    if (is.null(w)) w <- rep(1,length(all.data))
    #w <- w/w[1]

    for (i in 2:length(all.data)) {
#        pool$sequenced.base = find.max.of.2lists(pool$sequenced.base, 
#            all.data[[i]]$sequenced.base)
        pool$coverage <- pool$coverage + (w[i] * all.data[[i]]$coverage)
        pool$average.coverage <- pool$average.coverage + 
            (w[i] * all.data[[i]]$average.coverage)
        pool$base.with..10.coverage <- pool$base.with..10.coverage + 
            (w[i] * all.data[[i]]$base.with..10.coverage)
    }
    return(.removeChr(pool, remove.chrs))
}

.removeChr <- function(pool, remove.chrs=c()) {
    pool$chr <- gsub("^chrchr", "chr", pool$chr)
    idx <- pool$chr %in% remove.chrs
    pool$coverage[idx] <- NA
    pool$average.coverage[idx] <- NA
    pool$base.with..10.coverage[idx] <- NA
    pool
}
