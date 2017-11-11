context("correctCoverageBias")

normal.coverage.file <- system.file("extdata", "example_normal.txt", 
    package = "PureCN")
interval.file <- system.file("extdata", "ex2_intervals.txt", 
    package = "PureCN", mustWork = TRUE)
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
    package = "PureCN")

test_that("Example data matches after normalization", {
    output.file <- tempfile(fileext = ".txt")
    coverage <- correctCoverageBias(normal.coverage.file, gc.gene.file, 
        output.file = output.file)
    expect_equal(class(coverage)[1], "GRanges")
    expect_equal(length(coverage), 10049)
    correctCoverageBias(normal.coverage.file, gc.gene.file, plot.max.density = 100, 
        plot.gc.bias = TRUE)
    x <- readCoverageFile(output.file)
    expect_equal(x$average.coverage, coverage$average.coverage)
    correctCoverageBias(head(x, 200), gc.gene.file)
    gc.data <- read.delim(gc.gene.file, as.is = TRUE)
    gc.data$Gene <- NULL
    tmpFile <- tempfile()
    write.table(gc.data, file = tmpFile, row.names = FALSE, quote = FALSE, 
        sep = "\t")
    coverage2 <- correctCoverageBias(normal.coverage.file, tmpFile, 
        output.file = output.file)
    corCov <- cor(coverage$average.coverage, coverage2$average.coverage, 
        use = "complete.obs")
    expect_true(corCov > 0.99)
    output.png <- tempfile(fileext = ".png")
    png(file = output.png, width = 960)
    coverage <- correctCoverageBias(normal.coverage.file, gc.gene.file, 
        output.file = output.file, method = "POLYNOMIAL", 
        plot.gc.bias = TRUE)
    dev.off()
    file.remove(output.png)

    x <- readCoverageFile(output.file)
    expect_equal(x$average.coverage, coverage$average.coverage)
    file.remove(output.file)
})

test_that("Exceptions happen with wrong input", {
    expect_error(correctCoverageBias(normal.coverage.file, interval.file))
    expect_error(correctCoverageBias(normal.coverage.file, gc.gene.file, 
        method = "HELLOWORLD"))
})
