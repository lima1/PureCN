context("readSegmentationFile")

test_that("Example DNAcopy data matches", {
    seg.file <- system.file("extdata", "example_seg.txt", 
        package = "PureCN")
    seg <- readSegmentationFile(seg.file, "Sample1")
    offset <- -0.0033
    expect_equal(54, nrow(seg))
    expect_equal(0.133381833060556 - offset, seg$seg.mean[1], tolerance = .0001)
    expect_equal(-0.6394 - offset, seg$seg.mean[54], tolerance = .0001)
})

test_that("Example GATK4 data matches", {
    seg.file <- system.file("extdata", "example_gatk4_modelfinal.seg.gz", 
        package = "PureCN")
    seg <- readSegmentationFile(seg.file, "Sample1")
    offset <- -0.0037
    expect_equal(23, nrow(seg))
    expect_equal(-0.004295 - offset, seg$seg.mean[1], tolerance = .0001)
    expect_equal(0.002534 - offset, seg$seg.mean[23], tolerance = .0001)
})
