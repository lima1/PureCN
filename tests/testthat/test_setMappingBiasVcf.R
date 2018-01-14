context("setMappingBiasVcf")

vcf.file <- system.file("extdata", "example_vcf.vcf", package = "PureCN")
vcf <- readVcf(vcf.file, "hg19")

test_that("Mapping bias without normal panel matches", {
    vcf.bias <- round(setMappingBiasVcf(vcf)$bias, digits = 3)
    expected <- rep(0.977, 2331)
    expect_equal(vcf.bias, expected)
    vcf <- readVcf(vcf.file, "hg19", param = ScanVcfParam(samples = "LIB-02240e4"))
    vcf.bias <- round(setMappingBiasVcf(vcf)$bias, digits = 3)
    expected <- rep(1, 2331)
    expect_equal(vcf.bias, expected)
})

test_that("Mapping bias with normal panel matches", {
    normal_panel <- system.file("extdata", "normalpanel.vcf.gz", 
        package = "PureCN")
    mb <- setMappingBiasVcf(vcf, normal.panel.vcf.file = normal_panel)
    idx <- mb$pon.count > 0
    expect_equal(head(mb$pon.count[idx], 3), c(15, 5, 27))
    expect_equal(head(mb$bias[idx], 3), c(0.3362525, 1.0002354, 
        1.0119481), tolerance = 0.001)
    nvcf <- readVcf(normal_panel, genome = "hg19")
    idx <- overlapsAny(vcf, nvcf)
    expect_equal(sum(!mb$pon.count[idx] > 0), 0)
    expect_equal(sum(mb$pon.count[!idx] > 0), 0)
    expect_error(setMappingBiasVcf(vcf, normal.panel.vcf.file = normal_panel, 
        min.normals = 1))
})

test_that("Precomputed mapping bias matches", {
    normal_panel <- system.file("extdata", "normalpanel.vcf.gz", 
        package = "PureCN")
    mb <- calculateMappingBiasVcf(normal_panel)
    ov <- findOverlaps(vcf, mb, select = "first")
    idx <- !is.na(ov)
    expect_equivalent(head(mb$pon.count[ov[idx]], 3), c(15, 5, 27))
    expect_equal(head(mb$bias[ov[idx]],3), c(0.3362525, 1.0002354, 
        1.0119481), tolerance = 0.001)

    normal_panel_precomp <- tempfile(fileext = ".rds")
    saveRDS(mb, file = normal_panel_precomp)
    mb <- setMappingBiasVcf(vcf, normal.panel.vcf.file = normal_panel_precomp)
    idx <- mb$pon.count > 0
    expect_equal(head(mb$pon.count[idx], 3), c(15, 5, 27))
    expect_equal(head(mb$bias[idx], 3), c(0.3362525, 1.0002354, 
        1.0119481), tolerance = 0.001)
    file.remove(normal_panel_precomp)
})
