#' Multi sample normalization and segmentation
#'
#' This function performs normalization and segmentation when multiple
#' for the same patient are available.
#'
#' CURRENTLY DEFUNCT BECAUSE IT DEPENDS ON THE DEFUNCT COPYNUMBER PACKAGE.
#' We are working on a replacement.
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
#' normal1.coverage.file <- system.file("extdata", "example_normal.txt.gz",
#'     package = "PureCN")
#' normal2.coverage.file <- system.file("extdata", "example_normal2.txt.gz",
#'     package = "PureCN")
#' tumor1.coverage.file <- system.file("extdata", "example_tumor.txt.gz",
#'     package = "PureCN")
#' tumor2.coverage.file <- system.file("extdata", "example_tumor2.txt.gz",
#'     package = "PureCN")
#'
#' normal.coverage.files <- c(normal1.coverage.file, normal2.coverage.file)
#' tumor.coverage.files <- c(tumor1.coverage.file, tumor2.coverage.file)
#'
#' normalDB <- createNormalDatabase(normal.coverage.files)
#'
#' # seg <- processMultipleSamples(tumor.coverage.files,
#' #         sampleids = c("Sample1", "Sample2"),
#' #         normalDB = normalDB,
#' #         genome = "hg19")
#'
#' @export processMultipleSamples
processMultipleSamples <- function(tumor.coverage.files, sampleids, normalDB,
    num.eigen = 20, genome, plot.cnv = TRUE, w = NULL,
    min.interval.weight = 1 / 3,
    max.segments = NULL, chr.hash = NULL, centromeres = NULL, ...) {
    .Defunct(msg="preprocessMultipleSamples is temporarily defunct")
}
