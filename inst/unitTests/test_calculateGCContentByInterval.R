test_calculateGCContentByInterval <- function() {
    reference.file <- system.file("extdata", "ex2_reference.fa", 
        package = "PureCN", mustWork = TRUE)
    interval.file <- system.file("extdata", "ex2_intervals.txt", 
        package = "PureCN", mustWork = TRUE)
    bed.file <- system.file("extdata", "ex2_intervals.bed", 
        package="PureCN", mustWork = TRUE)
    gc <- calculateGCContentByInterval(interval.file, reference.file, 
        output.file = "gc_file.txt")
    x <- read.delim("gc_file.txt")
    checkEqualsNumeric(c(0.4533333,0.5057143,0.5733333,0.4800000,0.3600000), x$gc_bias, tolerance=0.001) 
    checkEqualsNumeric(c(0.4533333,0.5057143,0.5733333,0.4800000,0.3600000), gc$gc_bias, tolerance=0.001) 
    intervals <- import(bed.file)
    y <- calculateGCContentByInterval(intervals, reference.file, 
        output.file="gc_file_bed_test.txt")
    checkEqualsNumeric(x$gc_bias, y$gc_bias)
    checkEqualsNumeric(x$Target, as.character(y))
    y <- read.delim("gc_file_bed_test.txt")
    checkEqualsNumeric(x$gc_bias, y$gc_bias)

    interval.file2 <- "ex2_intervals2.txt"
    idata <- read.delim(interval.file, as.is=TRUE)
    idata[3,1] <- "seq2:0-149"
    write.table(idata, file=interval.file2, row.names=FALSE, quote=FALSE)
    checkException(calculateGCContentByInterval(interval.file2, reference.file))
    checkTrue(grepl("Interval coordinates should start at 1, not at 0", geterrmessage()))
    gc <- calculateGCContentByInterval(interval.file, reference.file, off.target=TRUE, min.off.target.width=2, off.target.padding=-2)
    checkEqualsNumeric(11, length(gc))
}    
