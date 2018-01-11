context("bootstrapResults")

test_that("Bootstrapping removed solutions", {
    data(purecn.example.output)
    set.seed(123)
    ret <- bootstrapResults(purecn.example.output, n = 100, top = 2)
    expect_equal(ret$results[[1]]$purity, purecn.example.output$results[[1]]$purity)
    expect_equal(ret$results[[1]]$ploidy, purecn.example.output$results[[1]]$ploidy)
    expect_true(length(ret$results) < length(purecn.example.output$results))
    expect_true(ret$results[[1]]$bootstrap.value >= 0.5)
    expect_true(ret$results[[2]]$bootstrap.value < 0.5)
    expect_true(length(ret$results) >= 2)
    ret <- bootstrapResults(purecn.example.output, n = 100, top = 3)
    expect_true(length(ret$results) >= 3)
    ret <- bootstrapResults(purecn.example.output, n = 100)
    expect_equal(length(purecn.example.output$results), length(ret$results))
})

