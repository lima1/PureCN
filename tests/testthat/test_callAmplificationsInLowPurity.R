context("callAmplificationsInLowPurity")

normal.coverage.file <- system.file("extdata", "example_normal.txt", 
    package = "PureCN")
normal2.coverage.file <- system.file("extdata", "example_normal2.txt", 
    package = "PureCN")
normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)

test_that("Example is called correctly", {
    data(purecn.example.output)
    normalDB <- createNormalDatabase(normal.coverage.files)
    m <- callAmplificationsInLowPurity(purecn.example.output,
       normalDB, all.genes = TRUE)

    esr2 <- m["ESR2", ]
    expect_equal(as.character(esr2$chr), "chr14")
    expect_true(esr2$start > 64694600)
    expect_true(esr2$end < 64761128)
})
