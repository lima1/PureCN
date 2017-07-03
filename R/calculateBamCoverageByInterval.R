#' Function to calculate coverage from BAM file
#' 
#' Takes a BAM file and an interval file as input and returns coverage for each
#' interval. Coverage should be then GC-normalized using the
#' \code{\link{correctCoverageBias}} function before determining purity and
#' ploidy with \code{\link{runAbsoluteCN}}. Uses the \code{scanBam} function
#' and applies low quality, duplicate reads as well as secondary alignment
#' filters.
#' 
#' 
#' @param bam.file Filename of a BAM file.
#' @param interval.file File specifying the intervals. Interval is expected in
#' first column in format CHR:START-END. The \code{gc.gene.file} can be used.
#' @param output.file Optionally, write minimal coverage file. Can be read with
#' the \code{\link{readCoverageFile}} function.
#' @param index.file The bai index. This is expected without the .bai file
#' suffix, see \code{?scanBam}.
#' @param keep.duplicates Keep or remove duplicated reads.
#' @return Returns total and average coverage by intervals.
#' @author Markus Riester
#' @seealso \code{\link{calculateGCContentByInterval}
#' \link{correctCoverageBias} \link{runAbsoluteCN}}
#' @examples
#' 
#' bam.file <- system.file("extdata", "ex1.bam", package="PureCN", 
#'     mustWork = TRUE)
#' interval.file <- system.file("extdata", "ex1_intervals.txt", 
#'     package="PureCN", mustWork = TRUE)
#' 
#' # Calculate raw coverage from BAM file. These need to be corrected for GC-bias
#' # using the correctCoverageBias function before determining purity and ploidy.
#' coverage <- calculateBamCoverageByInterval(bam.file=bam.file, 
#'     interval.file=interval.file)
#' 
#' @export calculateBamCoverageByInterval
#' @importFrom Rsamtools headerTabix ScanBamParam scanBamFlag
#'             scanBam scanFa scanFaIndex TabixFile
calculateBamCoverageByInterval <- function(bam.file, interval.file, 
    output.file = NULL, index.file = bam.file, keep.duplicates = FALSE) {
    intervalGr <- readCoverageFile(interval.file)

    param <- ScanBamParam(what=c("pos", "qwidth"), 
                which=intervalGr, 
                flag=scanBamFlag(isUnmappedQuery=FALSE, 
                             isNotPassingQualityControls=FALSE, 
                             isSecondaryAlignment=FALSE,
                             isDuplicate=if (keep.duplicates) NA else FALSE
                     )
             )

    x <- scanBam(bam.file, index=index.file, param=param)
    intervalGr$coverage <- sapply(seq_along(x), function(i) 
        sum(coverage(IRanges(x[[i]][["pos"]], width=x[[i]][["qwidth"]]), 
            shift=-start(intervalGr)[i], width=width(intervalGr)[i] )))

    intervalGr$average.coverage <- 
        intervalGr$coverage/width(intervalGr)

    if (!is.null(output.file)) {
        .writeCoverage(intervalGr, output.file)
    }    
    invisible(intervalGr)
}

.writeCoverage <- function(intervalGr, output.file) {
    tmp <- data.frame(
        Target=as.character(intervalGr), 
        total_coverage=intervalGr$coverage, 
        average_coverage=intervalGr$average.coverage,
        on_target=intervalGr$on.target
    )
    write.table(tmp, file=output.file, row.names=FALSE, quote=FALSE)
}
