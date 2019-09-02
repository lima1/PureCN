#' Read file containing interval-level log2 tumor/normal ratios 
#' 
#' Read log2 ratio file produced by external tools like The Genome Analysis 
#' Toolkit version 4.
#' 
#' @param file Log2 coverage file.
#' @param format File format. If missing, derived from the file 
#' extension. Currently GATK4 DenoiseReadCounts format supported.
#' A simple GATK3-style format, two columns with coordinates
#' as string in format chr:start-stop in first and log2-ratio
#' in second is also supported.
#' @param zero Start position is 0-based. Default is \code{FALSE}
#' for GATK, \code{TRUE} for BED file based intervals.
#' @return A \code{GRange} with the log2 ratio.
#' @author Markus Riester
#' @examples
#' 
#' logratio.file <- system.file("extdata", "example_gatk4_denoised_cr.tsv.gz",
#'     package = "PureCN")
#' logratio <- readLogRatioFile(logratio.file)
#' 
#' @export readLogRatioFile
readLogRatioFile <- function(file, format, zero = NULL) {
    if (missing(format)) format <- .getLogRatioFormat(file)
    if (format == "GATK3") return(.readLogRatioFileGATK3(file, zero = FALSE))
    if (format == "GATK4") return(.readLogRatioFileGATK4(file, zero = FALSE))
}

.getLogRatioFormat <- function(file) {
     header <- scan(file, what = character(), sep = "\n", nmax = 1, quiet = TRUE)
     format <- "GATK4"
     if (grepl("^Target", header)[1]) return("GATK3")
     format
}

.readLogRatioFileGATK3 <- function(file, zero=FALSE) {
    x <- fread(file, data.table = FALSE)
    gr <- GRanges(x[,1])
    gr$log.ratio <- x[,2]
    gr
}

.readLogRatioFileGATK4 <- function(file, zero = FALSE) {
    x <- read.delim(file, comment.char = "@", as.is = TRUE)
    gr <- GRanges(x[,1], IRanges(start = x[,2], end = x[,3]))
    gr$log.ratio <- x[,4]
    gr
}
