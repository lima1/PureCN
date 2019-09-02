context("readLogRatioFile")

test_that("Example data matches", {
    logratio.file <- system.file("extdata", "example_gatk4_denoised_cr.tsv.gz", 
        package = "PureCN")
    logratio <- readLogRatioFile(logratio.file)
    expect_equal(21, length(logratio))
    expect_equal(0.109473, logratio$log.ratio[1], tolerance = .00001)
    expect_equal(-0.185664, logratio$log.ratio[21], tolerance = .00001)
    logratio.file2 <- system.file("extdata", "example_logratio.txt.gz", 
        package = "PureCN")
    logratio2 <- readLogRatioFile(logratio.file2)
    expect_equal(as.character(logratio), as.character(logratio2))
    expect_equal(logratio$log.ratio, logratio2$log.ratio)
})
