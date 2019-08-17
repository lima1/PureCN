context("readLogRatioFile")

test_that("Example data matches", {
    logratio.file <- system.file("extdata", "example_gatk4_denoised_cr.tsv.gz", 
        package = "PureCN")
    logratio <- readLogRatioFile(logratio.file)
    expect_equal(21, length(logratio))
    expect_equal(0.109473, logratio$log.ratio[1], tolerance = .00001)
    expect_equal(-0.185664, logratio$log.ratio[21], tolerance = .00001)
})
