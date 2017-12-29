context("calculateBamCoverageByInterval")

output.file <- tempfile(fileext = ".txt")

test_that("Coverage from test BAM file matches", {
    bam.file <- system.file("extdata", "ex1.bam", package = "PureCN", 
        mustWork = TRUE)
    interval.file <- system.file("extdata", "ex1_intervals.txt", 
        package = "PureCN", mustWork = TRUE)
    coverage <- calculateBamCoverageByInterval(bam.file = bam.file, 
        interval.file = interval.file, output.file = output.file)
    expect_equal(coverage$average.coverage, c(20.95205, 43.78357, 
        21.29271), tolerance=0.01)
    expect_equal(coverage$counts, c(610, 1158, 636), tolerance=0.01)
    expect_equal(unlist(coverage$duplication.rate), rep(0, 3), check.names=FALSE)
}) 

test_that("Coverage output is correct", {
    x <- readCoverageFile(output.file)
    expect_equal(x$average.coverage, c(20.95205, 43.78357, 21.29271), tolerance=0.01)
    expect_equal(x$counts, c(610, 1158, 636), tolerance=0.01)
    interval.file <- system.file("extdata", "example_intervals.txt", 
        package = "PureCN")
    expect_error(correctCoverageBias(x, interval.file))
})

file.remove(output.file)
