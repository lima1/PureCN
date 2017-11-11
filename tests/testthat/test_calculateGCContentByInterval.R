test_that("test_calculateGCContentByInterval", {
    reference.file <- system.file("extdata", "ex2_reference.fa", 
        package = "PureCN", mustWork = TRUE)
    interval.file <- system.file("extdata", "ex2_intervals.txt", 
        package = "PureCN", mustWork = TRUE)
    bed.file <- system.file("extdata", "ex2_intervals.bed", package = "PureCN", 
        mustWork = TRUE)
    gc <- calculateGCContentByInterval(interval.file, reference.file, 
        output.file = "gc_file.txt")
    x <- read.delim("gc_file.txt", as.is=TRUE)
    expect_equal(x$gc_bias, c(0.4533333, 0.5057143, 0.5733333, 
        0.48, 0.36), tolerance=0.001)
    expect_equal(gc$gc_bias, c(0.4533333, 0.5057143, 0.5733333, 
        0.48, 0.36), tolerance=0.001)
    intervals <- import(bed.file)
    y <- calculateGCContentByInterval(intervals, reference.file, 
        output.file = "gc_file_bed_test.txt")
    expect_equal(y$gc_bias, x$gc_bias)
    expect_equal(as.character(y), x$Target)
    y <- read.delim("gc_file_bed_test.txt")
    expect_equal(y$gc_bias, x$gc_bias)
    interval.file2 <- "ex2_intervals2.txt"
    idata <- read.delim(interval.file, as.is = TRUE)
    idata[3, 1] <- "seq2:0-149"
    write.table(idata, file = interval.file2, row.names = FALSE, 
        quote = FALSE)
    expect_error(calculateGCContentByInterval(interval.file2, 
        reference.file))
    expect_true(grepl("Interval coordinates should start at 1, not at 0", 
        geterrmessage()))
    gc <- calculateGCContentByInterval(interval.file, reference.file, 
        off.target = TRUE, min.off.target.width = 2, off.target.padding = -2)
    expect_equal(length(gc), 11)
    intervals <- import(bed.file)
    gc2 <- calculateGCContentByInterval(gc, reference.file)
    expect_equal(start(gc2), start(gc))
    expect_equal(end(gc2), end(gc))
    expect_equal(gc2$mappability, gc$mappability)
    expect_equal(gc2$gc_bias, gc$gc_bias)
    if (.Platform$OS.type != "windows") {
        mappability.file <- system.file("extdata", "ex2_mappability.bigWig", 
            package = "PureCN", mustWork = TRUE)
        mappability <- import(mappability.file)
        gcMap <- calculateGCContentByInterval(intervals, reference.file, 
            mappability = mappability)
        expect_equal(gcMap$mappability, c(1, 1, 0.7, 1, 1), tolerance = 0.001)
    }
    mappability.file <- system.file("extdata", "ex2_mappability.bed", 
        package = "PureCN", mustWork = TRUE)
    mappability <- import(mappability.file)
    gcMap <- calculateGCContentByInterval(intervals, reference.file, 
        mappability = mappability)
    expect_equal(gcMap$mappability, c(1, 1, 0.7, 1, 1), tolerance = 0.001)
    reference.file <- system.file("extdata", "ex3_reference.fa", 
        package = "PureCN", mustWork = TRUE)
    bed.file3 <- system.file("extdata", "ex3_intervals.bed", 
        package = "PureCN", mustWork = TRUE)
    intervals3 <- import(bed.file3)
    x <- calculateGCContentByInterval(intervals3, reference.file)
    expect_equal(x$gc_bias, c(0.4533333, 0.5057143, 0.5733333,
        0.48, 0.36), tolerance = 0.001)
    seqlevelsStyle(intervals3) <- "NCBI"
    x <- calculateGCContentByInterval(intervals3, reference.file)
    expect_equal(x$gc_bias, c(0.4533333, 0.5057143, 0.5733333, 
        0.48, 0.36), tolerance = 0.001)
    expect_error(calculateGCContentByInterval(intervals, reference.file))
    expect_true(grepl("Chromosome naming style of interval file", 
        geterrmessage()))
    mappability.file3 <- system.file("extdata", "ex3_mappability.bed", 
        package = "PureCN", mustWork = TRUE)
    mappability3 <- import(mappability.file3)
    seqlevelsStyle(mappability3) <- "NCBI"
    x <- calculateGCContentByInterval(intervals3, reference.file, 
        mappability = mappability3)
    expect_equal(x$gc_bias, c(0.4533333, 0.5057143, 0.5733333, 
        0.48, 0.36), tolerance = 0.001)
    expect_equal(x$mappability, c(1, 1, 0.7, 1, 1), tolerance = 0.001)
})
