context("setMappingBiasVcf")

vcf.file <- system.file("extdata", "example_vcf.vcf.gz", package = "PureCN")
vcf <- readVcf(vcf.file, "hg19")

test_that("Mapping bias without normal panel matches", {
    vcf.bias <- round(info(setMappingBiasVcf(vcf))$MBB, digits = 3)
    expected <- rep(0.977, 2331)
    expect_equal(vcf.bias, expected)
    vcf <- readVcf(vcf.file, "hg19", param = ScanVcfParam(samples = "LIB-02240e4"))
    vcf.bias <- round(info(setMappingBiasVcf(vcf))$MBB, digits = 3)
    expected <- rep(1, 2331)
    expect_equal(vcf.bias, expected)
})

test_that("Precomputed mapping bias matches", {
    normal_panel <- system.file("extdata", "normalpanel.vcf.gz", 
        package = "PureCN")
    mb <- calculateMappingBiasVcf(normal_panel, genome = "hg19")
    ov <- findOverlaps(vcf, mb, select = "first")
    idx <- !is.na(ov)
    expect_equivalent(head(mb$pon.count[ov[idx]], 3), c(15, 5, 27))
    expect_equal(head(mb$bias[ov[idx]],3), c(0.3362525, 1.0002354, 
        1.0119481), tolerance = 0.001)

    normal_panel_precomp <- tempfile(fileext = ".rds")
    saveRDS(mb, file = normal_panel_precomp)
    mb <- setMappingBiasVcf(vcf, mapping.bias.file = normal_panel_precomp)
    idx <- info(mb)$MBPON > 0
    expect_equal(head(info(mb)$MBPON[idx], 3), c(15, 5, 27))
    expect_equal(head(info(mb)$MBB[idx], 3), c(0.3362525, 1.0002354, 
        1.0119481), tolerance = 0.001)
    expect_equal( weighted.mean(c(1.00023541351692, 1.01194812157128, 0.915037402634247), c(5,27,12)),
        info(mb)$MBB[230], tolerance = 0.001)
    file.remove(normal_panel_precomp)

    vcf.single.file <- system.file("extdata", "example_single.vcf.gz", package = "PureCN")
    expect_error(calculateMappingBiasVcf(vcf.single.file), "only a single sample")
    expect_error(setMappingBiasVcf(vcf, mapping.bias.file = normal_panel), "rds suffix")
})

test_that("GenomicsDB import works", {
    skip_if_not(requireNamespace("genomicsdb"), "genomicsdb required")
    skip_if_not(requireNamespace("jsonlite"), "jsonlite required")
    resources_file <- system.file("extdata", "gatk4_pon_db.tgz", 
        package = "PureCN")
    tmp_dir <- tempdir()
    untar(resources_file, exdir = tmp_dir)
    workspace <- file.path(tmp_dir, "gatk4_pon_db")
    bias <- calculateMappingBiasGatk4(workspace, "hg19")
    expect_equal(2101, length(bias))
    unlink(tmp_dir, recursive=TRUE)
})


test_that("Issue 184_2 is fixed", {
    vcf.184.2 <- readVcf(system.file("extdata", "issue184_2.vcf.gz",
        package = "PureCN"))
    mb.184.2 <- readRDS(system.file("extdata", "issue184_2_mb.rds",
        package = "PureCN"))
    expect_equal(1, PureCN:::.findOverlapsCheckAlt(vcf.184.2, mb.184.2))
    expect_equal(2, PureCN:::.findOverlapsCheckAlt(vcf.184.2, rev(mb.184.2)))
    alt(vcf.184.2) <- DNAStringSetList("A")
    expect_equal(2, PureCN:::.findOverlapsCheckAlt(vcf.184.2, mb.184.2))
    expect_equal(1, PureCN:::.findOverlapsCheckAlt(vcf.184.2, rev(mb.184.2)))
})    
