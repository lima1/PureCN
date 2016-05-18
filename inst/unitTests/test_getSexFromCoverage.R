test_getSexFromCoverage <- function() {
    gatk.tumor.file <- system.file("extdata", "example_tumor.txt", package="PureCN")
    coverage <- readCoverageGatk(gatk.tumor.file)
    sex <- getSexFromCoverage(coverage, verbose=FALSE)
    checkTrue(is.na(sex))
    sex <- getSexFromCoverage(gatk.tumor.file, verbose=FALSE)
    checkTrue(is.na(sex))

    chr22 <- coverage[which(coverage$chr=="chr22"),]
    chrX <- chr22
    chrX$chr <- "chrX"
    chrY <- chr22
    chrY$chr <- "chrY"
    coverage.sex <- do.call(rbind, list(coverage, chrX, chrY))
    sex <- getSexFromCoverage(coverage.sex)
    checkIdentical(sex, "M")

    chrY$average.coverage <- chrY$average.coverage/50
    coverage.sex <- do.call(rbind, list(coverage, chrX, chrY))
    sex <- getSexFromCoverage(coverage.sex)
    checkIdentical(sex, "F")

    chrY <- chr22
    chrY$chr <- "chrY"
    chrY$average.coverage <- chrY$average.coverage/10
    coverage.sex <- do.call(rbind, list(coverage, chrX, chrY))
    sex <- getSexFromCoverage(coverage.sex)
    checkTrue(is.na(sex))
}    
