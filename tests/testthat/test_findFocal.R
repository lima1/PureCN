context("findFocal")

test_that("Example data is called correctly", {
    data(purecn.example.output)
    ret <- findFocal(purecn.example.output$results[[1]]$seg)
    expect_equal(class(ret), "logical")
    expect_true(nrow(purecn.example.output$results[[1]]$seg) == 
        length(ret))
    expect_true(min(purecn.example.output$results[[1]]$seg[ret, 
        "C"]) >= 5)
})
