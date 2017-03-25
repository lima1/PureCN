test_readCoverageFile <- function() {
    tumor.coverage.file <- system.file("extdata", "example_tumor.txt", package="PureCN")
    coverage <- readCoverageFile(tumor.coverage.file)
    checkEquals(nrow(coverage), 10049)
    checkIdentical(as.character(coverage$chr[1]), "chr1")

    tumor.overlapping.coverage.file <- system.file("extdata", "test_coverage_overlapping_intervals.txt", package="PureCN")
    coverage <- readCoverageFile(tumor.overlapping.coverage.file)
    checkEquals(5, nrow(coverage))
    checkEqualsNumeric(c(1216042, 1216045, 1216606, 1216791, 1216991), coverage$probe_start)
}    
