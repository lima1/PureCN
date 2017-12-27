context("calculateGCContentByInterval")

reference.file <- system.file("extdata", "ex2_reference.fa", 
    package = "PureCN", mustWork = TRUE)
interval.file <- system.file("extdata", "ex2_intervals.txt", 
    package = "PureCN", mustWork = TRUE)
bed.file <- system.file("extdata", "ex2_intervals.bed", package = "PureCN", 
    mustWork = TRUE)

output.file <- tempfile(fileext = ".txt")

test_that("GC-bias of example reference and intervals (GATK format) matches", {
    gc <- calculateGCContentByInterval(interval.file, reference.file, 
        output.file = output.file)
    x <- read.delim(output.file, as.is = TRUE)
    expect_equal(x$gc_bias, c(0.4533333, 0.5057143, 0.5733333, 
        0.48, 0.36), tolerance = 0.001)
    expect_equal(gc$gc_bias, c(0.4533333, 0.5057143, 0.5733333, 
        0.48, 0.36), tolerance = 0.001)
})

test_that("GC-bias of example reference and intervals (BED format) matches", {
    x <- read.delim(output.file, as.is = TRUE)
    intervals <- import(bed.file)
    output.file2 <- tempfile(fileext = ".txt")
    y <- calculateGCContentByInterval(intervals, reference.file, 
        output.file = output.file2)
    expect_equal(y$gc_bias, x$gc_bias)
    expect_equal(as.character(y), x$Target)
    y <- read.delim(output.file2)
    file.remove(output.file2)
    expect_equal(y$gc_bias, x$gc_bias)
})

test_that("Exceptions happen with wrong input", {
    interval.file2 <- tempfile(fileext = ".txt")
    idata <- read.delim(interval.file, as.is = TRUE)
    idata[3, 1] <- "seq2:0-149"
    write.table(idata, file = interval.file2, row.names = FALSE, 
        quote = FALSE)
    expect_error(calculateGCContentByInterval(interval.file2, 
        reference.file),
        "Interval coordinates should start at 1, not at 0")
    file.remove(interval.file2)
})

test_that("reptiming annotated correctly", {
    reptiming.file <- system.file("extdata", "ex2_reptiming.bed",
        package = "PureCN", mustWork = TRUE)

    reptiming <- import(reptiming.file)
    intervals <- import(bed.file)
    gr <- calculateGCContentByInterval(intervals, reference.file, 
        reptiming = reptiming)
    expect_equal(c(17.5, 11.0, 50, 10.0, 10.0), gr$reptiming)
})

test_that("Offtarget settings work as expected", {
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
    expect_error(calculateGCContentByInterval(intervals, reference.file),
        "Chromosome naming style of interval file")
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

file.remove(output.file)
