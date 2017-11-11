test_that("test_readCoverageFile", {
    tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
        package = "PureCN")
    coverage <- readCoverageFile(tumor.coverage.file)
    expect_equal(length(coverage), 10049)
    expect_identical("chr1", as.character(seqnames(coverage)[1]))
    expect_equal(sum(!is.na(coverage[seqnames(coverage) == "chr21"]$coverage)), 
        179)
    pool <- poolCoverage(list(coverage))
    expect_equal(pool$average.coverage, coverage$average.coverage)
    pool <- poolCoverage(list(coverage), remove.chrs = "chr21")
    expect_equal(sum(is.na(pool[seqnames(coverage) == "chr21"]$coverage)), 
        179)
    tumor.overlapping.coverage.file <- system.file("extdata", 
        "test_coverage_overlapping_intervals.txt", package = "PureCN")
    coverage <- readCoverageFile(tumor.overlapping.coverage.file)
    expect_equal(length(coverage), 3)
    expect_equal(start(coverage), c(1216042, 1216606, 1216791))
    expect_equal(end(coverage), c(1216050, 1216678, 1217991))
    coverageFile <- system.file("extdata", "example_normal3.cnn", 
        package = "PureCN")
    coverage <- readCoverageFile(coverageFile)
    expect_equal(length(coverage), 4)
    expect_equal(start(coverage), c(762097, 861281, 865591, 866325) + 
        1)
    expect_equal(end(coverage), c(762270, 861490, 865791, 866498))
    expect_equal(coverage$on.target, c(TRUE, TRUE, TRUE, TRUE))
    coverage <- readCoverageFile(coverageFile, zero = FALSE)
    expect_equal(length(coverage), 4)
    expect_equal(start(coverage), c(762097, 861281, 865591, 866325))
    expect_equal(end(coverage), c(762270, 861490, 865791, 866498))
    expect_equal(coverage$on.target, c(TRUE, TRUE, TRUE, TRUE))
    coverageFile <- system.file("extdata", "example_normal4.cnr", 
        package = "PureCN")
    coverage <- readCoverageFile(coverageFile)
    expect_equal(length(coverage), 5)
    expect_equal(start(coverage), c(10500, 70509, 227917, 318219, 
        367658) + 1)
    expect_equal(end(coverage), c(68590, 176917, 267219, 367158, 
        367893))
    expect_equal(coverage$on.target, c(FALSE, FALSE, FALSE, FALSE, 
        TRUE))
})

