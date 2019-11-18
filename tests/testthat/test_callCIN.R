context("callCIN")

test_that("Example is called correctly", {
    data(purecn.example.output)
    loh <- callLOH(purecn.example.output)
    loh$size <- loh$end - loh$start + 1
    idx <- loh$C == 2
    ret <- callCIN(purecn.example.output, reference.state = "normal",
                   allele.specific = FALSE)
    expect_equal(sum(loh$size[!idx])/sum(loh$size), ret, tol = 0.001)
    loh <- loh[!is.na(loh$M),]
    ret <- callCIN(purecn.example.output)
    expect_equal(0.481, ret, tol = 0.01)
    ret <- callCIN(purecn.example.output, reference.state = "normal")
    idx <- loh$C == 2 & loh$M == 1
    expect_equal(sum(loh$size[!idx])/sum(loh$size), ret, tol = 0.001)
})
