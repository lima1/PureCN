test_findFocal <- function() {
    gatk.normal.file <- system.file("extdata", "example_normal.txt", 
        package="PureCN")
    gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
        package="PureCN")
    coverage <- correctCoverageBias(gatk.normal.file, gc.gene.file, 
        output.file="test_loess_coverage.txt")

    checkEquals("data.frame", class(coverage))
    checkEquals(10049, nrow(coverage))

    checkTrue(nrow(purecn.example.output$results[[1]]$seg)==length(ret))
    checkTrue( min(purecn.example.output$results[[1]]$seg[ret,"C"]) >= 6)

    x <- readCoverageGatk("test_loess_coverage.txt")
    checkEquals(c(20.95205,43.78357,21.29271), x$average.coverage, tolerance=0.01)
    checkEqualsNumeric(coverage$average.coverage, x$average.coverage)
    file.remove("test_loess_coverage.txt")
}
