context("callAlterations")

test_that("Example is called correctly", {
    data(purecn.example.output)
    calls <- callAlterations(purecn.example.output)
    expect_true(sum(calls$C < 6 & calls$C > 0.5) == 0)
    calls <- callAlterations(purecn.example.output, failed = TRUE)
    expect_true(sum(calls$gene.mean < 0.9 & calls$gene.mean >
        -0.9) == 0)
    esr2 <- callAlterations(purecn.example.output, all.genes = TRUE)["ESR2", ]
    expect_equal(as.character(esr2$chr), "chr14")
    expect_true(esr2$start > 64694600)
    expect_true(esr2$end < 64761128)
})

test_that("issue_292 is fixed", {
    data(purecn.example.output)
    calls <- callAlterations(purecn.example.output, id = 2, all.genes = TRUE)
    expect_true(abs(mean(calls$C) - purecn.example.output$results[[2]]$ploidy) < 0.5)
})
