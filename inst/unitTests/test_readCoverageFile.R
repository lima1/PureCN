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

    coverageFile <- system.file("extdata", "example_normal3.cnn", package="PureCN")
    coverage <- readCoverageFile(coverageFile)
    checkEquals(4, length(coverage))
    checkEqualsNumeric(c(762097, 861281, 865591, 866325)+1, start(coverage))
    checkEqualsNumeric(c(762270, 861490, 865791, 866498), end(coverage))
    checkEquals(c(TRUE, TRUE, TRUE, TRUE), coverage$on.target)

    coverage <- readCoverageFile(coverageFile, zero=FALSE)
    checkEquals(4, length(coverage))
    checkEqualsNumeric(c(762097, 861281, 865591, 866325), start(coverage))
    checkEqualsNumeric(c(762270, 861490, 865791, 866498), end(coverage))
    checkEquals(c(TRUE, TRUE, TRUE, TRUE), coverage$on.target)

    coverageFile <- system.file("extdata", "example_normal4.cnr", package="PureCN")
    coverage <- readCoverageFile(coverageFile)
    checkEquals(5, length(coverage))
    checkEqualsNumeric(c(10500, 70509, 227917, 318219, 367658)+1, start(coverage))
    checkEqualsNumeric(c(68590, 176917, 267219, 367158, 367893), end(coverage))
    checkEquals(c(FALSE, FALSE, FALSE, FALSE, TRUE), coverage$on.target)
}    
