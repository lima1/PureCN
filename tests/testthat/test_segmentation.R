context("segmentation")

normal.coverage.file <- system.file("extdata", "example_normal_tiny.txt",
    package = "PureCN")
tumor.coverage.file <- system.file("extdata", "example_tumor_tiny.txt",
    package = "PureCN")
vcf.file <- system.file("extdata", "example.vcf.gz",
    package = "PureCN")
seg.file <- system.file("extdata", "example_seg.txt",
    package = "PureCN")

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
   skip_if_not(PureCN:::.checkGATK4Version("4.1.7.0") >= 0,
    "gatk binary > 4.1.7.0 required")

    ret <- runAbsoluteCN(normal.coverage.file = normal.coverage.file,
        tumor.coverage.file = tumor.coverage.file, vcf.file = vcf.file,
        sampleid = "Sample1",  genome = "hg19",
        fun.segmentation = segmentationGATK4, max.ploidy = 4,
        test.purity = seq(0.3, 0.7, by = 0.05),
        max.candidate.solutions = 1, plot.cnv = FALSE)

    expect_equal(0.65, ret$results[[1]]$purity, tolerance = 0.02)
    expect_equal(1.62, ret$results[[1]]$ploidy, tolerance = 0.2)
})

test_that("Hclust segmentation works", {
    expect_error(runAbsoluteCN(normal.coverage.file = normal.coverage.file,
        tumor.coverage.file = tumor.coverage.file,
        sampleid = "Sample1",  genome = "hg19",
        fun.segmentation = segmentationHclust,
        max.candidate.solutions = 1, plot.cnv = FALSE),
        "segmentationHclust requires an")
})


test_that("private function .fixBreakpoint.", {
    seg <- readSegmentationFile(seg.file, "Sample1")
    data(purecn.example.output)
    gr <- purecn.example.output$input$log.ratio
    lr <- gr$log.ratio
    seg_1 <- PureCN:::.fixBreakpointsInBaits(gr, lr, seg, purecn.example.output$input$chr.hash)
    expect_equivalent(seg_1$loc.start, seg$loc.start)
    expect_equivalent(seg_1$loc.end, seg$loc.end)

    seg[24, "loc.start"] <- 82403793 + 1
    seg[44, "loc.end"] <- 57507347

    seg_1 <- PureCN:::.fixBreakpointsInBaits(gr, lr, seg, purecn.example.output$input$chr.hash)

    expect_equivalent(seg[23, "loc.start"], seg_1[23, "loc.start"])
    expect_equivalent(82403838, seg_1[23, "loc.end"])
    expect_equivalent(82403838 + 1, seg_1[24, "loc.start"])
    expect_equivalent(seg[24, "loc.end"], seg_1[24, "loc.end"])

    expect_equivalent(seg[44, "loc.start"], seg_1[44, "loc.start"])
    expect_equivalent(57507289 - 1, seg_1[44, "loc.end"])
    expect_equivalent(57507289, seg_1[45, "loc.start"])
    expect_equivalent(seg[45, "loc.end"], seg_1[45, "loc.end"])

    expect_equivalent(seg$loc.start[-c(23, 24, 44, 45)],
        seg_1$loc.start[-c(23, 24, 44, 45)])
    expect_equivalent(seg$loc.end[-c(23, 24, 44, 45)],
        seg_1$loc.end[-c(23, 24, 44, 45)])
})

test_that("issue 201 is fixed.", {
    expect_error(runAbsoluteCN(normal.coverage.file = normal.coverage.file,
        tumor.coverage.file = tumor.coverage.file,
        sampleid = "Sample1",  genome = "hg19",
        args.segmentation = list(undo.SD = "A"),
        max.candidate.solutions = 1, plot.cnv = FALSE),
        "undo.SD")
})    
