context("readLogRatioFile")
data(purecn.example.output)

test_that("Example data matches", {
    logratio.file <- system.file("extdata", "example_gatk4_denoised_cr.tsv.gz", 
        package = "PureCN")
    logratio <- readLogRatioFile(logratio.file)
    expect_equal(21, length(logratio))
    expect_equal(0.109473, logratio$log.ratio[1], tolerance = .00001)
    expect_equal(-0.185664, logratio$log.ratio[21], tolerance = .00001)
    expect_equivalent(seqlengths(logratio), c(248956422, 242193529, 156040895))
    logratio.file2 <- system.file("extdata", "example_logratio.txt.gz", 
        package = "PureCN")
    logratio2 <- readLogRatioFile(logratio.file2)
    expect_equal(as.character(logratio), as.character(logratio2))
    expect_equal(logratio$log.ratio, logratio2$log.ratio)
})

test_that("parsing -> writing -> parsing works", {
    x <- purecn.example.output$input
    y <- x
    y$log.ratio$log.ratio <- NULL
    output.file <- tempfile(fileext = ".tsv")
    expect_error(
        PureCN:::.writeLogRatioFileGATK4(y, 1, output.file),
        "log.ratio NULL"      
    )
    PureCN:::.writeLogRatioFileGATK4(x, 1, output.file)
    z <- readLogRatioFile(output.file)
    expect_equivalent(x$log.ratio$log.ratio, z$log.ratio)
    file.remove(output.file)
})

