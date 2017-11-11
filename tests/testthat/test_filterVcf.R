context("filterVcf")

vcf.file <- system.file("extdata", "example_vcf.vcf", package = "PureCN")
vcf <- readVcf(vcf.file, "hg19")

test_that("snp.blacklist filtering works", {
    x <- filterVcfBasic(vcf)
    d.f <- data.frame(id = head(names(vcf)), Count = 1)
    output.file <- tempfile(fileext = ".bed")
    rtracklayer::export(head(rowRanges(vcf)), con = output.file, 
        format = "bed")
    z <- filterVcfBasic(vcf, snp.blacklist = output.file)
    expect_equal(nrow(x$vcf) - nrow(z$vcf), 6)
    expect_equal(sum(d.f$id %in% names(x$vcf)), 6)
    expect_equal(sum(d.f$id %in% names(z$vcf)), 0)
    file.remove(output.file)
})

test_that("stats.file filtering works", {
    x <- filterVcfBasic(vcf)
    gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
        package = "PureCN")
    vcfMutectFilter <- filterVcfMuTect(vcf, stats.file = gc.gene.file)
    expect_equal(nrow(vcfMutectFilter$vcf), nrow(x$vcf))
    output.file <- tempfile(fileext = ".txt")
    cat("#testfile\n", file = output.file)
    d.f <- data.frame(head(rowRanges(vcf)))[, 1:2]
    colnames(d.f) <- c("contig", "position")
    suppressWarnings(write.table(d.f, file = output.file, append = TRUE, 
        quote = FALSE, row.names = FALSE, sep = "\t"))
    vcfMutectFilter <- filterVcfMuTect(vcf, stats.file = output.file)
    expect_equal(nrow(vcfMutectFilter$vcf), 6)
    file.remove(output.file)
})
