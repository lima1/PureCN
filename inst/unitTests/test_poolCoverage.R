test_poolCoverage <- function() {
    normal.coverage.file <- system.file("extdata", "example_normal.txt", 
        package = "PureCN")
    normal2.coverage.file <- system.file("extdata", "example_normal2.txt", 
        package = "PureCN")
    normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)
    coverage <- lapply(normal.coverage.files, readCoverageGatk)
    pool <- poolCoverage(coverage)
    checkEqualsNumeric( pool$average.coverage , coverage[[1]]$average.coverage+coverage[[2]]$average.coverage)
    checkEqualsNumeric( pool$coverage , coverage[[1]]$coverage+coverage[[2]]$coverage)

}    
