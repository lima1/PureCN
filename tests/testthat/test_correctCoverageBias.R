test_that("test_correctCoverageBias", {
    normal.coverage.file <- system.file("extdata", "example_normal.txt", 
        package = "PureCN")
    gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
        package = "PureCN")
    coverage <- correctCoverageBias(normal.coverage.file, gc.gene.file, 
        output.file = "test_loess_coverage.txt")
    expect_equal(class(coverage)[1], "GRanges")
    expect_equal(length(coverage), 10049)
    correctCoverageBias(normal.coverage.file, gc.gene.file, plot.max.density = 100, 
        plot.gc.bias = TRUE)
    x <- readCoverageFile("test_loess_coverage.txt")
    expect_equal(x$average.coverage, coverage$average.coverage)
    file.remove("test_loess_coverage.txt")
    interval.file <- system.file("extdata", "ex2_intervals.txt", 
        package = "PureCN", mustWork = TRUE)
    expect_error(correctCoverageBias(normal.coverage.file, interval.file))
    correctCoverageBias(head(x, 200), gc.gene.file)
    gc.data <- read.delim(gc.gene.file, as.is = TRUE)
    gc.data$Gene <- NULL
    tmpFile <- tempfile()
    write.table(gc.data, file = tmpFile, row.names = FALSE, quote = FALSE, 
        sep = "\t")
    coverage2 <- correctCoverageBias(normal.coverage.file, tmpFile, 
        output.file = "test_loess_coverage.txt")
    corCov <- cor(coverage$average.coverage, coverage2$average.coverage, 
        use = "complete.obs")
    expect_true(corCov > 0.99)
    png(file = "test_gc_bias.png", width = 960)
    coverage <- correctCoverageBias(normal.coverage.file, gc.gene.file, 
        output.file = "test_norm_coverage.txt", method = "POLYNOMIAL", 
        plot.gc.bias = TRUE)
    dev.off()
    file.remove("test_gc_bias.png")
    x <- readCoverageFile("test_norm_coverage.txt")
    expect_equal(x$average.coverage, coverage$average.coverage)
    file.remove("test_norm_coverage.txt")
    expect_error(correctCoverageBias(normal.coverage.file, gc.gene.file, 
        method = "HELLOWORLD"))
})

