test_readCoverageGatk <- function() {
    gatk.tumor.file <- system.file("extdata", "example_tumor.txt", package="PureCN")
    coverage <- readCoverageGatk(gatk.tumor.file)
    checkEquals(nrow(coverage), 10049)
    checkIdentical(as.character(coverage$chr[1]), "chr1")
}    
