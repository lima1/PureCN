context("segmentation")

test_that("Precomputed boudaries are correct", {
    data(purecn.DNAcopy.bdry)
    alpha <- formals(segmentationCBS)$alpha
    eta <- formals(segment)$eta
    nperm <- formals(segment)$nperm
    max.ones <- floor(nperm * alpha) + 1
    set.seed(123)
    sbdry <- getbdry(eta, nperm, max.ones)
    expect_equal(purecn.DNAcopy.bdry, sbdry)
})
