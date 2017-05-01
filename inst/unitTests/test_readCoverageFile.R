test_readCoverageFile <- function() {
    tumor.coverage.file <- system.file("extdata", "example_tumor.txt", package="PureCN")
    coverage <- readCoverageFile(tumor.coverage.file)
    checkEquals(10049, length(coverage))
    checkIdentical(as.character(seqnames(coverage)[1]), "chr1")

    tumor.overlapping.coverage.file <- system.file("extdata", "test_coverage_overlapping_intervals.txt", package="PureCN")
    coverage <- readCoverageFile(tumor.overlapping.coverage.file)
    checkEquals(3, length(coverage))
    checkEqualsNumeric(c(1216042, 1216606, 1216791), start(coverage))
    checkEqualsNumeric(c(1216050, 1216678, 1217991), end(coverage))
}    
