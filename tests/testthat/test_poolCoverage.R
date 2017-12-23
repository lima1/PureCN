context("poolCoverage")

normal.coverage.file <- system.file("extdata", "example_normal.txt", 
    package = "PureCN")
normal2.coverage.file <- system.file("extdata", "example_normal2.txt", 
    package = "PureCN")
normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)

test_that("Example coverage is averaged", {
    coverage <- lapply(normal.coverage.files, readCoverageFile)
    pool <- poolCoverage(coverage)
    expect_equal(coverage[[1]]$average.coverage + coverage[[2]]$average.coverage, 
        pool$average.coverage)
    expect_equal(coverage[[1]]$coverage + coverage[[2]]$coverage, 
        pool$coverage)
    pool2 <- poolCoverage(coverage, w=c(0.5, 0.5))
    expect_equal((coverage[[1]]$coverage+coverage[[2]]$coverage)/2, 
        pool2$coverage)
})

test_that("Exceptions happend with wrong input", {
    coverage <- lapply(normal.coverage.files, readCoverageFile)
    expect_error(poolCoverage(coverage, w=1:3), "different lengths")
})    
