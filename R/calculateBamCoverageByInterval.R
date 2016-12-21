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
#' the \code{\link{readCoverageGatk}} function.
#' @param index.file The bai index. This is expected without the .bai file
#' suffix, see \code{?scanBam}.
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
#'             scanBam scanFa TabixFile
calculateBamCoverageByInterval <- function(bam.file, interval.file, 
    output.file = NULL, index.file = bam.file) {
    interval <- read.delim(interval.file, as.is=TRUE)
    colnames(interval)[1] <- "Target"
    pos <- as.data.frame(do.call(rbind, strsplit(interval$Target, ":|-")), 
        stringsAsFactors = FALSE)
    interval.gr <- GRanges(seqnames = pos[,1], 
        IRanges(start = as.numeric(pos[,2]), end = as.numeric(pos[,3])))

    param <- ScanBamParam(what=c("pos", "qwidth"), 
                which=interval.gr, 
                flag=scanBamFlag(isUnmappedQuery=FALSE, 
                                 isNotPassingQualityControls=FALSE, 
                                 isSecondaryAlignment=FALSE,
                                 isDuplicate=FALSE))

    x <- scanBam(bam.file, index=index.file, param=param)
    cvg <- sapply(seq_along(x), function(i) 
        sum(coverage(IRanges(x[[i]][["pos"]], width=x[[i]][["qwidth"]]), 
            shift=-start(interval.gr)[i], width=width(interval.gr)[i] )))

    ret <- data.frame(
        probe=interval$Target, 
        chr=seqnames(interval.gr), 
        probe_start=start(interval.gr),
        probe_end=end(interval.gr),
        targeted.base=width(interval.gr),
        sequenced.base=NA,
        coverage=cvg
    )
    ret$average.coverage <- ret$coverage/ret$targeted.base
    ret$base.with..10.coverage <- NA

    if (!is.null(output.file)) {
        tmp <- ret[, c("probe", "coverage", "average.coverage")]
        colnames(tmp) <- c("Target", "total_coverage", "average_coverage")
        write.table(tmp, file=output.file, row.names=FALSE, quote=FALSE)
    }    
    invisible(ret)
}
