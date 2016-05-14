test_calculateGCContentByInterval <- function() {
    reference.file <- system.file("extdata", "ex2_reference.fa", 
        package = "PureCN", mustWork = TRUE)
    interval.file <- system.file("extdata", "ex2_intervals.txt", 
        package = "PureCN", mustWork = TRUE)
    gc <- calculateGCContentByInterval(interval.file, reference.file, 
        output.file = "gc_file.txt")
    x <- read.delim("gc_file.txt")
    checkEqualsNumeric(c(0.4533333,0.5057143,0.5733333,0.4800000,0.3600000), x$gc_bias, tolerance=0.001) 
    checkEqualsNumeric(c(0.4533333,0.5057143,0.5733333,0.4800000,0.3600000), gc$gc_bias, tolerance=0.001) 
}    

