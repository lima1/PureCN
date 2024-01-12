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

test_that("skipping base quality works", {
    f1 <- filterVcfBasic(vcf, min.base.quality = NULL)
    f2 <- filterVcfBasic(vcf, min.base.quality = 0)
    expect_equal(length(f1$vcf), length(f2$vcf))
    expect_error(filterVcfBasic(vcf, min.base.quality = 50), 
                 "No variants passed")
})

test_that("M2 VCF with POP_AF flag is annotated with DB flag", {
    # first check that the POP_AF field is parsed
    vcf.m2.file <- system.file("extdata", "example_mutect2.vcf.gz", package = "PureCN")
    vcf.m2 <- PureCN:::.readAndCheckVcf(vcf.m2.file, "hg38")
    expect_equal(c(TRUE, rep(FALSE, 10)), info(vcf.m2)$DB)
    expect_equal(unlist(info(vcf.m2)$POP_AF > 0.001), info(vcf.m2)$DB)

    expect_output(filterVcfMuTect(vcf.m2, use.somatic.status = FALSE), 
        "Less than half of variants are annoted as germline database member")
    output.file <- tempfile(fileext = ".vcf")
    writeVcf(vcf.m2, file = output.file)
    expect_output(vcf.m2 <- PureCN:::.readAndCheckVcf(output.file, "hg38"),
        "Will ignore POP_AF")
    expect_equal(c(TRUE, rep(FALSE, 10)), info(vcf.m2)$DB)
    expect_equal(unlist(info(vcf.m2)$POP_AF > 0.001), info(vcf.m2)$DB)
    file.remove(output.file)

    vcf.m2 <- PureCN:::.readAndCheckVcf(vcf.m2.file, "hg38")
    # testing dbSNP annotation only (rs* variant ids)
    # first create new VCF without the DB and POP_AF fields
    info(vcf.m2)$POP_AF <- NULL
    info(vcf.m2)$DB <- NULL
    info(header(vcf.m2)) <-  info(header(vcf.m2))[-match(c("DB", "POP_AF"), 
        rownames(info(header(vcf.m2)))),]
    output.file <- tempfile(fileext = ".vcf")
    writeVcf(vcf.m2, file = output.file)
    # now parse
    vcf.m2 <- PureCN:::.readAndCheckVcf(output.file, "hg38")
    file.remove(output.file)
    expect_equal(c(TRUE, rep(FALSE, 10)), info(vcf.m2)$DB)
    expect_equal(grepl("rs", names(vcf.m2)), info(vcf.m2)$DB)
})


test_that("issue 62 is fixed", {
    vcf.file <- system.file("extdata", "issue62.vcf.gz", package = "PureCN")
    x <- PureCN:::.readAndCheckVcf(vcf.file)
    expect_equivalent(c(901,53), as.vector(table(info(x)$SOMATIC)))
})              

test_that("issue 109 is fixed", {
    vcf.file <- system.file("extdata", "issue109.vcf.gz", package = "PureCN")
    expect_output(x <- PureCN:::.readAndCheckVcf(vcf.file), "AD field misses ref counts")
    expect_equivalent(c(272, 2), geno(x)$AD[[1,1]])
    expect_equivalent(274, geno(x)$DP[1,1])
    expect_equivalent(2/274, geno(x)$FA[[1,1]])
})

test_that("Missing FA does not cause crash", {
     tmp <- geno(vcf)$FA[2,1][[1]] 
     geno(vcf)$FA[2,1][[1]] <- NA
     expect_output( 
        x <- PureCN:::.readAndCheckVcf(vcf, "hg19"), rownames(vcf)[2])
     expect_equivalent(PureCN:::.countVariants(vcf),
                       PureCN:::.countVariants(x) + 1)
     geno(vcf)$FA[2,1][[1]] <- tmp
})    

test_that("issue 184 is fixed", {
    vcf.file <- system.file("extdata", "issue184.vcf.gz", package = "PureCN")
    vcf.184 <- PureCN:::.readAndCheckVcf(vcf.file)
    expect_output(x <- PureCN:::.getTumorIdInVcf(vcf.184), "GT field in VCF contains missing values")
    expect_equivalent("TC_098_1.2", x)
})

test_that("Missing BQ is handled correctly", {
    vcf.file <- system.file("extdata", "example_vcf.vcf.gz", package = "PureCN")
    vcf <- PureCN:::.readAndCheckVcf(vcf.file, "hg19")
    expect_equal(range(geno(vcf)$BQ[, 1]), c(8, 36))
    vcf.bq <- vcf
    geno(vcf.bq)$BQ <- NULL
    tmp <- rownames(geno(header(vcf.bq)))
    geno(header(vcf.bq)) <-  geno(header(vcf.bq))[-match("BQ", tmp), ]
    output.file <- tempfile(fileext = ".vcf")
    writeVcf(vcf.bq, file = output.file)
    vcf.bq2 <- PureCN:::.readAndCheckVcf(output.file, "hg19")
    expect_equal(range(geno(vcf.bq2)$BQ[, 1]), c(30, 30))
    file.remove(output.file)
})

test_that("Providing min.supporting.reads works", {
    x <- filterVcfBasic(vcf, use.somatic.status = FALSE, min.supporting.reads = 10)
    expect_equal(min(sapply(geno(x$vcf)$AD[,1], function(x) x[2])), 10)
})

test_that("issue 249", {
    # set random BQ to NA
    vcf.bq <- vcf
    geno(vcf.bq)$BQ[[2, 1]] <- NA
    vcf.bq <- PureCN:::.readAndCheckVcf(vcf.bq, "hg19")
    expect_equal(NA, geno(vcf.bq)$BQ[[2, 1]])
})

test_that("issue 320", {
    # set random BQ to NA
    expect_output(x <- filterVcfBasic(vcf, use.somatic.status = FALSE, min.base.quality = 34),
        "Many variants removed by min.base.quality")
    expect_equal(34, min(sapply(geno(x$vcf)$BQ[,1], function(x) x[1])))
})
