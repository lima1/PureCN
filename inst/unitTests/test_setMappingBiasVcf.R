test_setMappingBiasVcf <- function() {
    vcf.file <- system.file("extdata", "example_vcf.vcf", package = "PureCN")
    vcf <- readVcf(vcf.file, "hg19")
    vcf.bias <- round(setMappingBiasVcf(vcf)$bias, digits=3)
    expected <- rep(0.977, 2331)
    checkEquals(expected, vcf.bias)
    vcf <- readVcf(vcf.file, "hg19", param=ScanVcfParam(samples="LIB-02240e4"))
    vcf.bias <- round(setMappingBiasVcf(vcf)$bias, digits=3)
    expected <- rep(1, 2331)
    checkEquals(expected, vcf.bias)
    normal_panel <- system.file("extdata", "normalpanel.vcf.gz", package = "PureCN")
    vcf <- readVcf(vcf.file, "hg19")
    mb <- setMappingBiasVcf(vcf, normal.panel.vcf.file=normal_panel)
    idx <- mb$pon.count>0
    checkEqualsNumeric(c(15,5,27), head(mb$pon.count[idx],3))
    checkEqualsNumeric(c(0.3362525, 1.0002354, 1.0119481), head(mb$bias[idx],3), tolerance=0.01)
    nvcf <- readVcf(normal_panel, genome="hg19")
    idx <- overlapsAny(vcf, nvcf)
    # all overlapping have pon.count>0
    checkEqualsNumeric(0, sum(!mb$pon.count[idx]>0))
    # all non-overlapping have pon.count=0
    checkEqualsNumeric(0, sum(mb$pon.count[!idx]>0))

    checkException(setMappingBiasVcf(vcf, normal.panel.vcf.file=normal_panel, min.normals=1))
}    
