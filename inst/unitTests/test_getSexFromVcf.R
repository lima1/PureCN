test_getSexFromVcf <- function() {
    # a pretty useless test since the example data contains no homozygous calls.
    vcf.file <- system.file("extdata", "example_vcf.vcf", package = "PureCN")
    vcf <- readVcf(vcf.file, "hg19")
    sex <- getSexFromVcf(vcf)
    checkTrue(is.na(sex))
}    
