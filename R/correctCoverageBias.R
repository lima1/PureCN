correctCoverageBias <- structure(function(# Correct for GC bias
### Takes as input coverage data in GATK format (or data 
### read by \code{\link{readCoverageGatk}}) and a mapping file for 
### GC content, and then uses a loess normalization for bias correction. 
### Largely follows the GC correction of the TitanCNA package.
coverage.file, 
### Exon coverage file as produced by GATK. Either a file name
### or data parsed with the \code{\link{readCoverageGatk}} function.
gc.gene.file,
### File providing GC content for each exon in the coverage files.
### First column in format CHR:START-END. Second column GC content (0 to 1). 
### Third column provides gene symbols, which are optional, but used in
### \code{\link{runAbsoluteCN}} to generate gene level calls. This file
### can be generated with GATK GCContentByInterval tool or with the
### \code{\link{calculateGCContentByInterval}} function.
##seealso<< \code{\link{calculateGCContentByInterval}}
output.file=NULL
### Optionally, write file with GC corrected coverage. Can be read with
### the \code{\link{readCoverageGatk}} function.
) {
    if (is.character(coverage.file)) {
        tumor  <- readCoverageGatk(coverage.file)
    } else {
        tumor <- coverage.file
    }    
    
    gc <- .loadGcGeneFile(gc.gene.file, tumor)

    # taken from TitanCNA
    gc$valid <- TRUE
    gc$valid[tumor$average.coverage <= 0 | gc$gc_bias < 0] <- FALSE
    gc$ideal <- TRUE
    routlier <- 0.01
    range <- quantile(tumor$average.coverage[gc$valid], prob = 
        c(0, 1 - routlier), na.rm = TRUE)
    doutlier <- 0.001
    domain <- quantile(gc$gc_bias[gc$valid], prob = c(doutlier, 1 - doutlier), 
        na.rm = TRUE)
    gc$ideal[!gc$valid | 
        tumor$average.coverage <= range[1] |
        tumor$average.coverage > range[2] | 
        gc$gc_bias < domain[1] | 
        gc$gc_bias > domain[2]] <- FALSE

    rough <- loess(tumor$average.coverage[gc$ideal] ~ gc$gc_bias[gc$ideal], 
        span = 0.03)
    i <- seq(0, 1, by = 0.001)
    final <- loess(predict(rough, i) ~ i, span = 0.3)
    cor.gc <- predict(final, gc$gc_bias)
    cor.gc.factor <- cor.gc/mean(tumor$average.coverage, na.rm=TRUE)

    tumor$average.coverage <- tumor$average.coverage / cor.gc.factor
    tumor$coverage <- tumor$coverage / cor.gc.factor
    if (!is.null(output.file)) {
        tmp <- tumor[, c("probe", "coverage", "average.coverage")]
        colnames(tmp) <- c("Target", "total_coverage", "average_coverage")
        write.table(tmp, file=output.file, row.names=FALSE, quote=FALSE)
    }    
    invisible(tumor)
### GC normalized coverage.
}, ex=function() {
normal.coverage.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
    package="PureCN")
coverage <- correctCoverageBias(normal.coverage.file, gc.gene.file)
})

.loadGcGeneFile <- function(gc.gene.file, tumor) {
    gc <- read.delim(gc.gene.file)

    if (is.null(gc$gc_bias)) {
        .stopUserError("gc.gene.file header invalid.")
    }
    
    if (!identical(as.character(gc[,1]), as.character(tumor[,1]))) {
        if (sum(!as.character(tumor[,1]) %in% as.character(gc[,1])) > 0) {
            .stopUserError(
            "Interval files in coverage.file and gc.gene.file different.\n",
            "Some intervals in coverage have no GC information.")
        }
        warning(
        "Interval files in coverage.file and gc.gene.file different.")
        gc <- gc[match(as.character(tumor[,1]), as.character(gc[,1])),]
    }
    gc
}
    
