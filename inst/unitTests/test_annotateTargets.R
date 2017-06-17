test_annotateTargets <- function() {
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(org.Hs.eg.db)

    normal.coverage.file <- system.file("extdata", "example_normal.txt",
        package="PureCN")
    x <- head(readCoverageFile(normal.coverage.file),100)
    x <- annotateTargets(x,TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db)
    checkEquals("KIF1B", x$Gene[67])
}    
