context("processMultipleSamples")

test_that("example output correct", {
	normal1.coverage.file <- system.file("extdata", "example_normal.txt",
		package="PureCN")
	normal2.coverage.file <- system.file("extdata", "example_normal2.txt",
		package="PureCN")
	tumor1.coverage.file <- system.file("extdata", "example_tumor.txt",
		package="PureCN")
	tumor2.coverage.file <- system.file("extdata", "example_tumor2.txt",
		package="PureCN")

	normal.coverage.files <- c(normal1.coverage.file, normal2.coverage.file)
	tumor.coverage.files <- c(tumor1.coverage.file, tumor2.coverage.file)
	normalDB <- createNormalDatabase(normal.coverage.files)
    interval.weight.file <- tempfile(fileext = ".txt")
	calculateIntervalWeights(normal.coverage.files, interval.weight.file)

	seg <- processMultipleSamples(tumor.coverage.files,
			 sampleids = c("Sample1", "Sample2"),
			 normalDB = normalDB,
			 interval.weight.file = interval.weight.file,
			 genome = "hg19")
    expect_equal(c("Sample1", "Sample2"), levels(seg[,1]))
    seg.file <- tempfile(fileext = ".seg")
    write.table(seg, seg.file, row.names = FALSE, sep = "\t")
    vcf.file <- system.file("extdata", "example.vcf.gz", package = "PureCN")
    ret <- runAbsoluteCN(tumor.coverage.file = tumor1.coverage.file, 
        seg.file = seg.file, vcf.file = vcf.file, max.candidate.solutions = 1, 
        fun.segmentation = segmentationHclust,
        genome = "hg19", min.ploidy = 1.5, max.ploidy = 2.1, 
        test.purity = seq(0.4, 0.7, by = 0.05), sampleid = "Sample1")
    expect_equal(0.65, ret$results[[1]]$purity)

    file.remove(interval.weight.file)
})
