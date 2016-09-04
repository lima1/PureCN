test_setPriorVcf <- function() {
    vcf.file <- system.file("extdata", "example_vcf.vcf", package = "PureCN")
    vcf <- readVcf(vcf.file, "hg19")
    vcf.priorsomatic <- setPriorVcf(vcf)
    expected <- c(2322,9)
    names(expected) <- c(0.00001,0.999)
    checkEquals(expected[1], sort(table(vcf.priorsomatic))[2])
    checkEquals(expected[2], sort(table(vcf.priorsomatic))[1])
}    
