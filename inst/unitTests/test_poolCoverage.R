test_poolCoverage <- function() {
    gatk.normal.file <- system.file("extdata", "example_normal.txt", 
        package = "PureCN")
    gatk.normal2.file <- system.file("extdata", "example_normal2.txt", 
        package = "PureCN")
    gatk.normal.files <- c(gatk.normal.file, gatk.normal2.file)
    coverage <- lapply(gatk.normal.files, readCoverageGatk)
    pool <- poolCoverage(coverage)
    checkEqualsNumeric( pool$average.coverage , coverage[[1]]$average.coverage+coverage[[2]]$average.coverage)
    checkEqualsNumeric( pool$coverage , coverage[[1]]$coverage+coverage[[2]]$coverage)

}    
