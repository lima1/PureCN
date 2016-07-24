calculateBamCoverageByInterval <- structure(
function(# Function to calculate coverage from BAM file
### Takes a BAM file and an interval file as input and 
### returns coverage for each interval. Coverage should be GC-normalized
### using the \code{\link{correctCoverageBias}} function before determining
### purity and ploidy with \code{\link{runAbsoluteCN}}.
##seealso<< \code{\link{calculateGCContentByInterval} 
## \link{correctCoverageBias} \link{runAbsoluteCN}}
bam.file, 
### Filename of a BAM file.
interval.file,
### File specifying the intervals. Interval is expected in 
### first column in format CHR:START-END. The \code{gc.gene.file} can be used.
output.file=NULL
### Optionally, write minimal coverage file. Can be read with the
### \code{\link{readCoverageGatk}} function.
) {
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

    x <- scanBam(bam.file, param=param)
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
### Returns total and average coverage by intervals.
}, ex=function() {
bam.file <- system.file("extdata", "ex1.bam", package="PureCN", 
    mustWork = TRUE)
interval.file <- system.file("extdata", "ex1_intervals.txt", 
    package="PureCN", mustWork = TRUE)

coverage <- calculateBamCoverageByInterval(bam.file=bam.file, 
    interval.file=interval.file)
})    
