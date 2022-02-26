context("annotateTargets") 
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
test_coverage <- readCoverageFile(system.file("extdata", "example_normal.txt.gz",
    package = "PureCN"))

test_that("KIF1B is correctly annotated with UCSC chromosome names", {
    x <- head(test_coverage, 100)
    x <- annotateTargets(x, TxDb.Hsapiens.UCSC.hg19.knownGene, 
        org.Hs.eg.db)
    expect_equal(x$Gene[67], "KIF1B")
})

test_that("KIF1B is correctly annotated with NCBI chromosome names", {
    x <- head(test_coverage, 100)
    seqlevelsStyle(x) <- "Ensembl"
    x <- annotateTargets(x, TxDb.Hsapiens.UCSC.hg19.knownGene, 
        org.Hs.eg.db)
    expect_equal(x$Gene[67], "KIF1B")
})
