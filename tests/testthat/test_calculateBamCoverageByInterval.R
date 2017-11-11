test_that("test_calculateBamCoverageByInterval", {
    bam.file <- system.file("extdata", "ex1.bam", package = "PureCN", 
        mustWork = TRUE)
    interval.file <- system.file("extdata", "ex1_intervals.txt", 
        package = "PureCN", mustWork = TRUE)
    coverage <- calculateBamCoverageByInterval(bam.file = bam.file, 
        interval.file = interval.file, output.file = "ex1_coverage.txt")
    expect_equal(coverage$average.coverage, c(20.95205, 43.78357, 
        21.29271), tolerance=0.01)
    expect_equal(coverage$counts, c(610, 1158, 636), tolerance=0.01)
    expect_equal(unlist(coverage$duplication.rate), rep(0, 3), check.names=FALSE)
    x <- readCoverageFile("ex1_coverage.txt")
    expect_equal(x$average.coverage, c(20.95205, 43.78357, 21.29271), tolerance=0.01)
    expect_equal(x$counts, c(610, 1158, 636), tolerance=0.01)
    gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
        package = "PureCN")
    expect_error(correctCoverageBias(x, gc.gene.file))
})
