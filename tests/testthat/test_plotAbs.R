context("plotAbs")

test_that("Exceptions happen with wrong input", {
    data(purecn.example.output)
    expect_error( plotAbs(purecn.example.output, id = "hello", "BAF"), 
        "No solution with id hello")
    expect_error( plotAbs(purecn.example.output, id = 100, "BAF"), 
        "No solution with id 100")
})
