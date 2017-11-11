context("setPriorVcf")

test_that("Example data matches expected values", {
    vcf.file <- system.file("extdata", "example_vcf.vcf", package = "PureCN")
    vcf <- readVcf(vcf.file, "hg19")
    vcf.priorsomatic <- setPriorVcf(vcf)
    expected <- c(2322, 9)
    names(expected) <- c(1e-04, 0.999)
    expect_equal(sort(table(vcf.priorsomatic))[2], expected[1])
    expect_equal(sort(table(vcf.priorsomatic))[1], expected[2])
})

