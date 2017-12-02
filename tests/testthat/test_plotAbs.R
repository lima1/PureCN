context("plotAbs")

data(purecn.example.output)
ret <- predictSomatic(purecn.example.output)

test_that("Exceptions happen with wrong input", {
    expect_error( plotAbs(ret, id="hello", "BAF"), "No solution with id hello")
    expect_error( plotAbs(ret, id=10, "BAF"), "No solution with id 10")
})

