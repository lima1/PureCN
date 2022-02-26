context("correctCoverageBias")

normal.coverage.file <- system.file("extdata", "example_normal.txt.gz", 
    package = "PureCN")
interval.file <- system.file("extdata", "ex2_intervals.txt", 
    package = "PureCN", mustWork = TRUE)
interval.file2 <- system.file("extdata", "example_intervals.txt", 
    package = "PureCN")

test_that("Example data matches after normalization", {
    output.file <- tempfile(fileext = ".txt")
    coverage <- correctCoverageBias(normal.coverage.file, interval.file2, 
        output.file = output.file)
    expect_equal(class(coverage)[1], "GRanges")
    expect_equal(length(coverage), 10049)
    correctCoverageBias(normal.coverage.file, interval.file2, plot.max.density = 100, 
        plot.bias = TRUE)
    x <- readCoverageFile(output.file)
    expect_equal(x$average.coverage, coverage$average.coverage)
    correctCoverageBias(head(x, 200), interval.file2)
    gc.data <- read.delim(interval.file2, as.is = TRUE)
    gc.data$Gene <- NULL
    tmpFile <- tempfile()
    write.table(gc.data, file = tmpFile, row.names = FALSE, quote = FALSE, 
        sep = "\t")
    coverage2 <- correctCoverageBias(normal.coverage.file, tmpFile, 
        output.file = output.file)
    corCov <- cor(coverage$average.coverage, coverage2$average.coverage, 
        use = "complete.obs")
    expect_true(corCov > 0.99)
})

test_that("Exceptions happen with wrong input", {
    expect_error(correctCoverageBias(normal.coverage.file, interval.file))
    coverage <- readCoverageFile(normal.coverage.file)
    coverage$average.coverage <- 0
    expect_error(correctCoverageBias(coverage, interval.file), "zero")
})

test_that("Example data qc matches", {
    output.qc.file <- tempfile(fileext = ".txt")
    coverage <- correctCoverageBias(normal.coverage.file, interval.file2, 
        output.qc.file = output.qc.file)
    x <- read.delim(output.qc.file, sep=" ")
    expect_equal(1, nrow(x))
    expect_equal(10, ncol(x))
    file.remove(output.qc.file)
})

test_that("Example data without reptiming works", {
    x <- read.delim(interval.file2)
    interval.file3 <- tempfile(fileext = ".txt")
    x$reptiming <- NULL
    write.table(x, file=interval.file3, row.names=FALSE, quote=FALSE, sep="\t")
    coverage <- correctCoverageBias(normal.coverage.file, interval.file3, 
                                   plot.bias=TRUE)
    expect_equal(nrow(x), length(coverage))
    file.remove(interval.file3)
})
