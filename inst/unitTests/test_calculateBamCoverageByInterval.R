test_calculateBamCoverageByInterval <- function() {
    bam.file <- system.file("extdata", "ex1.bam", package="PureCN", 
        mustWork = TRUE)
    interval.file <- system.file("extdata", "ex1_intervals.txt", 
        package="PureCN", mustWork = TRUE)

    coverage <- calculateBamCoverageByInterval(bam.file=bam.file, 
        interval.file=interval.file, output.file="ex1_coverage.txt")
    
    checkEquals(c(20.95205,43.78357,21.29271), coverage$average.coverage, tolerance=0.01)
    x <- readCoverageGatk("ex1_coverage.txt")
    checkEquals(c(20.95205,43.78357,21.29271), x$average.coverage, tolerance=0.01)
}    
