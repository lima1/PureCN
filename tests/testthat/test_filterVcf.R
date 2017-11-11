test_that("test_filterVcf", {
    vcf.file <- system.file("extdata", "example_vcf.vcf", package = "PureCN")
    vcf <- readVcf(vcf.file, "hg19")
    x <- filterVcfBasic(vcf)
    d.f <- data.frame(id = head(names(vcf)), Count = 1)
    rtracklayer::export(head(rowRanges(vcf)), con = "snpbl2.bed", 
        format = "bed")
    z <- filterVcfBasic(vcf, snp.blacklist = "snpbl2.bed")
    expect_equal(nrow(x$vcf) - nrow(z$vcf), 6)
    expect_equal(sum(d.f$id %in% names(x$vcf)), 6)
    expect_equal(sum(d.f$id %in% names(z$vcf)), 0)
    gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
        package = "PureCN")
    vcfMutectFilter <- filterVcfMuTect(vcf, stats.file = gc.gene.file)
    expect_equal(nrow(vcfMutectFilter$vcf), nrow(x$vcf))
    filename <- "statstest.txt"
    cat("#testfile\n", file = filename)
    d.f <- data.frame(head(rowRanges(vcf)))[, 1:2]
    colnames(d.f) <- c("contig", "position")
    suppressWarnings(write.table(d.f, file = filename, append = TRUE, 
        quote = FALSE, row.names = FALSE, sep = "\t"))
    vcfMutectFilter <- filterVcfMuTect(vcf, stats.file = filename)
    expect_equal(nrow(vcfMutectFilter$vcf), 6)
})

