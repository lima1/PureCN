context("callLOH")

test_that("Example is called correctly", {
    data(purecn.example.output)
    ret <- callLOH(purecn.example.output)
    expect_true(is(ret, "data.frame"))
    expect_equal(13, ncol(ret))
})

test_that("NCBI-style chromosome names work", {
    normal.coverage.file <- system.file("extdata", "example_normal.txt", 
        package = "PureCN")
    tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
        package = "PureCN")
    vcf.file <- system.file("extdata", "example.vcf.gz", package = "PureCN")
    vcf <- readVcf(vcf.file)
    normal <- readCoverageFile(normal.coverage.file)
    tumor <- readCoverageFile(tumor.coverage.file)
    seqlevelsStyle(vcf) <- "Ensembl"
    seqlevelsStyle(normal) <- "Ensembl"
    seqlevelsStyle(tumor) <- "Ensembl"
    ret <- runAbsoluteCN(normal.coverage.file = normal, tumor.coverage.file = tumor, 
        genome = "hg19", vcf.file = vcf, sampleid = "Sample1", 
        min.ploidy = 1.4, max.ploidy = 2.4, test.purity = seq(0.4, 
            0.7, by = 0.05), max.candidate.solutions = 1, plot = FALSE)
    loh <- callLOH(ret)
    expect_equal(unique(loh$chr), as.character(1:22))
    loh2 <- callLOH(ret, keep.no.snp.segments = FALSE)
    expect_true(nrow(loh) > nrow(loh2))
    idx <- !is.na(loh$M)
    expect_equal(loh$C[idx], loh2$C)
    expect_equal(is.na(loh$M), is.na(loh$type))
})

test_that("No crash without centromeres", {
    x <- purecn.example.output
    x$input$centromeres <- NULL
    loh <- callLOH(x)
    expect_equal(13, ncol(loh))
})
