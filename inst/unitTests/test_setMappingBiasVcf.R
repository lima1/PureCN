test_setMappingBiasVcf <- function() {
    vcf.file <- system.file("extdata", "example_vcf.vcf", package = "PureCN")
    vcf <- readVcf(vcf.file, "hg19")
    vcf.bias <- round(setMappingBiasVcf(vcf), digits=3)
    expected <- rep(0.977, 2331)
    checkEquals(expected, vcf.bias)
    vcf <- readVcf(vcf.file, "hg19", param=ScanVcfParam(samples="LIB-02240e4"))
    vcf.bias <- round(setMappingBiasVcf(vcf), digits=3)
    expected <- rep(1, 2331)
    checkEquals(expected, vcf.bias)
}    
