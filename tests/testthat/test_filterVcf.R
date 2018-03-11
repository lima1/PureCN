context("filterVcf")

vcf.file <- system.file("extdata", "example_vcf.vcf.gz", package = "PureCN")
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
    interval.file <- system.file("extdata", "example_intervals.txt", 
        package = "PureCN")
    vcfMutectFilter <- filterVcfMuTect(vcf, stats.file = interval.file)
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

test_that("M2 VCF with POP_AF flag is annotated with DB flag", {
    # first check that the POP_AF field is parsed
    vcf.m2.file <- system.file("extdata", "example_mutect2.vcf.gz", package = "PureCN")
    vcf.m2 <- PureCN:::.readAndCheckVcf(vcf.m2.file, "hg38")
    expect_equal(c(TRUE, rep(FALSE, 10)), info(vcf.m2)$DB)
    expect_equal(unlist(info(vcf.m2)$POP_AF>0.001), info(vcf.m2)$DB)

    # testing dbSNP annotation only (rs* variant ids)
    # first create new VCF without the DB and POP_AF fields
    info(vcf.m2)$POP_AF <- NULL
    info(vcf.m2)$DB <- NULL
    info(header(vcf.m2)) <-  info(header(vcf.m2))[-match(c("DB", "POP_AF"), rownames(info(header(vcf.m2)))),]
    output.file <- tempfile(fileext = ".vcf")
    writeVcf(vcf.m2, file = output.file)
    # now parse
    vcf.m2 <- PureCN:::.readAndCheckVcf(output.file, "hg38")
    file.remove(output.file)
    expect_equal(c(TRUE, rep(FALSE, 10)), info(vcf.m2)$DB)
    expect_equal(grepl("rs", names(vcf.m2)), info(vcf.m2)$DB)
})
