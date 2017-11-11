test_that("test_getSexFromCoverage", {
    tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
        package = "PureCN")
    coverage <- readCoverageFile(tumor.coverage.file)
    sex <- getSexFromCoverage(coverage)
    expect_true(is.na(sex))
    sex <- getSexFromCoverage(tumor.coverage.file)
    expect_true(is.na(sex))
    chr22 <- coverage[which(seqnames(coverage) == "chr22")]
    chrX <- renameSeqlevels(chr22, c(chr22 = "chrX"))
    chrY <- renameSeqlevels(chr22, c(chr22 = "chrY"))
    coverage.sex <- c(coverage, chrX, chrY)
    sex <- getSexFromCoverage(coverage.sex)
    expect_identical("M", sex)
    chrY$average.coverage <- chrY$average.coverage/50
    coverage.sex <- c(coverage, chrX, chrY)
    sex <- getSexFromCoverage(coverage.sex)
    expect_identical("F", sex)
    chrY <- chr22
    chrY <- renameSeqlevels(chr22, c(chr22 = "chrY"))
    chrY$average.coverage <- chrY$average.coverage/21
    coverage.sex <- c(coverage, chrX, chrY)
    sex <- getSexFromCoverage(coverage.sex)
    expect_true(is.na(sex))
})

