#' Read file containing interval-level log2 tumor/normal ratios 
#' 
#' Read log2 ratio file produced by external tools like The Genome Analysis 
#' Toolkit version 4.
#' 
#' @param file Log2 coverage file.
#' @param format File format. If missing, derived from the file 
#' extension. Currently GATK4 DenoiseReadCounts format supported.
#' @param zero Start position is 0-based. Default is \code{FALSE}
#' for GATK, \code{TRUE} for BED file based intervals.
#' @return A \code{GRange} with the log2 ratio.
#' @author Markus Riester
#' @examples
#' 
#' tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
#'     package="PureCN")
#' coverage <- readCoverageFile(tumor.coverage.file)
#' 
#' @export readLogRatioFile
readLogRatioFile <- function(file, format, zero=NULL) {
    if (missing(format)) format <- "GATK4"
    if (format == "GATK4") return(.readLogRatioFileGATK4(file, zero = FALSE))
}

.readLogRatioFileGATK4 <- function(file, zero = FALSE) {
    x <- read.delim(file, comment.char = "@", as.is = TRUE)
    gr <- GRanges(x[,1], IRanges(start = x[,2], end = x[,3]))
    gr$log.ratio <- x[,4]
    gr
}
