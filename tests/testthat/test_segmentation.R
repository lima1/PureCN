context("segmentation")

normal.coverage.file <- system.file("extdata", "example_normal_tiny.txt", 
    package="PureCN")
tumor.coverage.file <- system.file("extdata", "example_tumor_tiny.txt", 
    package="PureCN")
vcf.file <- system.file("extdata", "example.vcf.gz",
    package="PureCN")

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


test_that("GATK4 wrapper works for example data.", {
   skip_if_not(PureCN:::.checkGATK4Version("4.1.7.0") >= 0, "gatk binary > 4.1.7.0 required")
 
    ret <-runAbsoluteCN(normal.coverage.file = normal.coverage.file, 
        tumor.coverage.file = tumor.coverage.file, vcf.file = vcf.file, 
        sampleid = "Sample1",  genome = "hg19",
        fun.segmentation = segmentationGATK4, max.ploidy = 4,
        test.purity = seq(0.3, 0.7, by = 0.05),
        max.candidate.solutions = 1, plot.cnv = FALSE)
    
    expect_equal(0.65, ret$results[[1]]$purity, tolerance = 0.02)
    expect_equal(1.62, ret$results[[1]]$ploidy, tolerance = 0.2)
})
