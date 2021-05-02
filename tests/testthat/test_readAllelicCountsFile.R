context("readAllelicCountsFile")

vcf.file <- system.file("extdata", "example.vcf.gz", package = "PureCN")
ac.file <- system.file("extdata", "example_allelic_counts.tsv", package = "PureCN")
vcf <- readVcf(vcf.file, "hg19")
data(purecn.example.output)
normal.coverage.file <- system.file('extdata', 'example_normal.txt', 
     package='PureCN')
tumor.coverage.file <- system.file('extdata', 'example_tumor.txt', 
     package='PureCN')

test_that("example parses correctly", {
    vcf_ac <- readAllelicCountsFile(ac.file)
    expect_equal(as.character(ref(vcf_ac)), as.character(ref(head(vcf,20))))
})

test_that("parsing -> writing -> parsing works", {
    output.file <- tempfile(fileext = ".tsv")
    PureCN:::.writeAllelicCountsFileGatk(vcf, 1, output.file)
    vcf_ac <- readAllelicCountsFile(output.file)
    expect_equal(as.character(ref(vcf_ac)), as.character(ref(vcf)))
    ret <- runAbsoluteCN(normal.coverage.file = normal.coverage.file,
        tumor.coverage.file = tumor.coverage.file,
        candidates = purecn.example.output$candidates,
        vcf.file = vcf,
        genome = "hg19",
        test.purity = seq(0.4, 0.7, by = 0.05), min.ploidy = 1.5, 
        max.ploidy = 2.4, max.candidate.solutions = 1, plot.cnv=FALSE)
    expect_true(length(ret$results)>0)
    file.remove(output.file)
})

