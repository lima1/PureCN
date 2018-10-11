context("getSexFromVcf")

test_that("Example data is called correctly", {
    vcf.file <- system.file("extdata", "example.vcf.gz", package = "PureCN")
    vcf <- readVcf(vcf.file, "hg19")
    sex <- getSexFromVcf(vcf)
    expect_true(is.na(sex))
    vcfs <- vcf[info(vcf)$SOMATIC]
    getSexFromVcf(vcfs, "LIB-02240e4")
    expect_true(is.na(sex))
})
