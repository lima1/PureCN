test_readCoverageGatk <- function() {
    tumor.coverage.file <- system.file("extdata", "example_tumor.txt", package="PureCN")
    coverage <- readCoverageGatk(tumor.coverage.file)
    checkEquals(nrow(coverage), 10049)
    checkIdentical(as.character(coverage$chr[1]), "chr1")
}    
