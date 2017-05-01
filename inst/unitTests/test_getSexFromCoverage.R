test_getSexFromCoverage <- function() {
    tumor.coverage.file <- system.file("extdata", "example_tumor.txt", package="PureCN")
    coverage <- readCoverageFile(tumor.coverage.file)
    sex <- getSexFromCoverage(coverage)
    checkTrue(is.na(sex))
    sex <- getSexFromCoverage(tumor.coverage.file)
    checkTrue(is.na(sex))

    chr22 <- coverage[which(seqnames(coverage)=="chr22")]
    chrX <- renameSeqlevels(chr22, c(chr22="chrX"))
    chrY <- renameSeqlevels(chr22, c(chr22="chrY"))
    coverage.sex <- c(coverage, chrX, chrY)
    sex <- getSexFromCoverage(coverage.sex)
    checkIdentical(sex, "M")

    chrY$average.coverage <- chrY$average.coverage/50
    coverage.sex <- c(coverage, chrX, chrY)
    sex <- getSexFromCoverage(coverage.sex)
    checkIdentical(sex, "F")

    chrY <- chr22
    chrY <- renameSeqlevels(chr22, c(chr22="chrY"))
    chrY$average.coverage <- chrY$average.coverage/21
    coverage.sex <- c(coverage, chrX, chrY)
    sex <- getSexFromCoverage(coverage.sex)
    checkTrue(is.na(sex))
}    
