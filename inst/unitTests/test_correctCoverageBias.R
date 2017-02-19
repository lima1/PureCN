test_correctCoverageBias <- function() {
    normal.coverage.file <- system.file("extdata", "example_normal.txt", 
        package="PureCN")
    gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
        package="PureCN")
    coverage <- correctCoverageBias(normal.coverage.file, gc.gene.file, 
        output.file="test_loess_coverage.txt")

    checkEquals("data.frame", class(coverage))
    checkEquals(10049, nrow(coverage))

    x <- readCoverageGatk("test_loess_coverage.txt")
    checkEqualsNumeric(coverage$average.coverage, x$average.coverage)
    file.remove("test_loess_coverage.txt")
    interval.file <- system.file("extdata", "ex2_intervals.txt", 
        package = "PureCN", mustWork = TRUE)
    checkException(correctCoverageBias(normal.coverage.file, interval.file))

    data(purecn.example.output)
    coveragePass2 <- correctCoverageBias(normal.coverage.file, gc.gene.file, 
        output.file="test_loess_coverage.txt", 
        purecn.output=purecn.example.output)
    corPass2 <- cor(coverage$average.coverage, coveragePass2$average.coverage, use="complete.obs")
    checkTrue(corPass2 > 0.99)
    correctCoverageBias(head(x,200), gc.gene.file)

    png(file="test_gc_bias.png",width=960)
    coverage <- correctCoverageBias(normal.coverage.file, gc.gene.file, 
        output.file="test_norm_coverage.txt", method = "POLYNOMIAL", plot.gc.bias = TRUE)
    dev.off()
    file.remove("test_gc_bias.png")
    x <- readCoverageGatk("test_norm_coverage.txt")
    checkEqualsNumeric(coverage$average.coverage, x$average.coverage)
    file.remove("test_norm_coverage.txt")
}
