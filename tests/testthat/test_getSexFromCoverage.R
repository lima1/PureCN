context("getSexFromCoverage")

tumor.coverage.file <- system.file("extdata", "example_tumor.txt.gz", 
    package = "PureCN")
coverage <- readCoverageFile(tumor.coverage.file)
chr22 <- coverage[which(seqnames(coverage) == "chr22")]
chrX <- renameSeqlevels(chr22, c(chr22 = "chrX"))

test_that("Warning with missing coverage data", {
    sex <- getSexFromCoverage(coverage) 
    expect_true(is.na(sex))
    expect_output( getSexFromCoverage(coverage), "WARN" )
})

test_that("Warning with missing coverage data in file", {
    sex <- getSexFromCoverage(tumor.coverage.file) 
    expect_true(is.na(sex))
    expect_output( getSexFromCoverage(coverage), "WARN" )
})

test_that("Male correct from coverage data", {
    chrY <- renameSeqlevels(chr22, c(chr22 = "chrY"))
    coverage_fakemale <- suppressWarnings(c(coverage, chrX, chrY))
    sex <- getSexFromCoverage(coverage_fakemale)
    expect_identical("M", sex)
})

test_that("Female correct from coverage data", {
    chrY <- renameSeqlevels(chr22, c(chr22 = "chrY"))
    chrY$average.coverage <- chrY$average.coverage/50
    coverage_fakefemale <- suppressWarnings( c(coverage, chrX, chrY) )
    sex <- getSexFromCoverage(coverage_fakefemale)
    expect_identical("F", sex)
})
    
test_that("NA correct from contaminated coverage data", {
    chrY <- renameSeqlevels(chr22, c(chr22 = "chrY"))
    chrY$average.coverage <- chrY$average.coverage / 21
    coverage_fakecontamination <- suppressWarnings( c(coverage, chrX, chrY))
    sex <- getSexFromCoverage(coverage_fakecontamination)
    expect_true(is.na(sex))
})
