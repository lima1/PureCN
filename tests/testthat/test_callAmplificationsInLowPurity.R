context("callAmplificationsInLowPurity")

test_that("Example is called correctly", {
    data(purecn.example.output)
     m <- callAmplificationsInLowPurity(purecn.example.output,
        normalDB, all.genes = TRUE)

    esr2 <- m["ESR2", ]
    expect_equal(as.character(esr2$chr), "chr14")
    expect_true(esr2$start > 64694600)
    expect_true(esr2$end < 64761128)
})
