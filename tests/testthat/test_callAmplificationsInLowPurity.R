context("callAmplificationsInLowPurity")

normal.coverage.file <- system.file("extdata", "example_normal.txt", 
    package = "PureCN")
normal2.coverage.file <- system.file("extdata", "example_normal2.txt", 
    package = "PureCN")
normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)
normalDB <- createNormalDatabase(normal.coverage.files)
data(purecn.example.output)

test_that("Example is called correctly", {
    m <- callAmplificationsInLowPurity(purecn.example.output,
       normalDB, all.genes = TRUE, purity = 0.65)
    m2 <- callAmplificationsInLowPurity(purecn.example.output,
       normalDB, all.genes = TRUE, purity = 0.65, BPPARAM=BiocParallel::bpparam())
    esr2 <- m["ESR2", ]
    expect_equal(as.character(esr2$chr), "chr14")
    expect_true(esr2$start > 64694600)
    expect_true(esr2$end < 64761128)
    expect_true(esr2$C < 3 && esr2$C >= 2)
    expect_gt(cor(m$p.value, m2$p.value), 0.99)
})
test_that("Exceptions happen with incorrect input data", {
    expect_error(callAmplificationsInLowPurity(purecn.example.output, 
        normalDB, pvalue.cutoff = 1.2), "pvalue.cutoff")
    expect_error(callAmplificationsInLowPurity(purecn.example.output, 
        normalDB, pvalue.cutoff = -1.2), "pvalue.cutoff")
    expect_error(callAmplificationsInLowPurity(purecn.example.output, 
        normalDB, percentile.cutoff = 120), "percentile.cutoff")
    expect_error(callAmplificationsInLowPurity(purecn.example.output, 
        normalDB, percentile.cutoff = -120), "percentile.cutoff")
    expect_error(callAmplificationsInLowPurity(purecn.example.output, 
        normalDB, purity = -120), "purity")
    expect_error(callAmplificationsInLowPurity(purecn.example.output, 
        normalDB, purity = 80), "purity")
})
