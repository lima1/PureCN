#' Find best normal sample in database
#' 
#' Function to find the best matching normal for a provided tumor sample.
#' This function is deprecated and most features are not relevant with
#' the replacement function \code{\link{calculateTangentNormal}}.
#' 
#' 
#' @param tumor.coverage.file Coverage file or data of a tumor sample.
#' @param normalDB Database of normal samples, created with
#' \code{\link{createNormalDatabase}}.
#' @param pcs Principal components to use for distance calculation.
#' @param num.normals Return the \code{num.normals} best normals.
#' @param ignore.sex If \code{FALSE}, detects sex of sample and returns best
#' normals with matching sex.
#' @param sex Sex of sample. If \code{NULL}, determine with
#' \code{\link{getSexFromCoverage}} and default parameters. Valid values are
#' \code{F} for female, \code{M} for male. If all chromosomes are diploid,
#' specify \code{diploid}.
#' @param normal.coverage.files Only consider these normal samples. If
#' \code{NULL}, use all in the database. Must match
#' \code{normalDB$normal.coverage.files}.
#' @param pool If \code{TRUE}, use \code{\link{poolCoverage}} to pool best
#' normals.
#' @param pool.weights Either find good pooling weights by optimization or
#' weight all best normals equally.
#' @param plot.pool Allows the pooling function to create plots.
#' @param \dots Additional arguments passed to \code{\link{poolCoverage}}.
#' @return Filename of the best matching normal.
#' @author Markus Riester
#' @seealso \code{\link{createNormalDatabase} \link{getSexFromCoverage}}
#' @export findBestNormal
#' @importFrom stats dist predict
findBestNormal <- function(tumor.coverage.file, normalDB, pcs=1:3, 
    num.normals = 1, ignore.sex = FALSE, sex = NULL, 
    normal.coverage.files = NULL, pool = FALSE, 
    pool.weights = c("voom", "equal"), plot.pool = FALSE, 
    ...) {
    .Deprecated("calculateTangentNormal")
    calculateTangentNormal(tumor.coverage.file, normalDB, ignore.sex = ignore.sex, 
        sex = sex)
}


#' Plot the PCA of tumor and its best normal(s)
#' 
#' This method is defunct with no replacement
#' 
#' @export plotBestNormal
plotBestNormal <- function() {
    .Defunct()
}
