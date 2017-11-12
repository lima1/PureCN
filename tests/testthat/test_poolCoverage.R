context("poolCoverage")

test_that("Example coverage is averaged", {
    normal.coverage.file <- system.file("extdata", "example_normal.txt", 
        package = "PureCN")
    normal2.coverage.file <- system.file("extdata", "example_normal2.txt", 
        package = "PureCN")
    normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)
    coverage <- lapply(normal.coverage.files, readCoverageFile)
    pool <- poolCoverage(coverage)
    expect_equal(coverage[[1]]$average.coverage + coverage[[2]]$average.coverage, 
        pool$average.coverage)
    expect_equal(coverage[[1]]$coverage + coverage[[2]]$coverage, 
        pool$coverage)
})
