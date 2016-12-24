library('getopt')

### Parsing command line ------------------------------------------------------

spec <- matrix(c(
'help' , 'h', 0, "logical",
'force' , 'f', 0, "logical",
'bam', 'b', 1, "character",
'gatkcoverage', 'g', 1, "character",
'gcgene' , 'c', 1, "character",
'method' , 'm', 1, "character",
'outdir' , 'o', 1, "character"
), byrow=TRUE, ncol=4)
opt <- getopt(spec)

if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE))
    q(status=1)
}

force <- !is.null(opt$force)

bam.file <- opt$bam
gatk.coverage <- opt$gatkcoverage
gc.gene.file <- opt$gcgene
outdir <- opt$outdir

method <- ifelse(is.null(opt$method), "LOESS", opt$method)

outdir <- normalizePath(outdir, mustWork=TRUE)
gc.gene.file <- normalizePath(gc.gene.file, mustWork=TRUE)

### Calculate coverage from BAM files -----------------------------------------

if (!is.null(bam.file)) {
    library(PureCN)
    bam.file <- normalizePath(bam.file, mustWork=TRUE)
    output.file <- file.path(outdir,  gsub(".bam$","_coverage.txt", 
        basename(bam.file)))

    if (file.exists(output.file) && !force) {
        message(output.file, " exists. Skipping... (--force will overwrite)")
    } else {
        calculateBamCoverageByInterval(bam.file=bam.file, 
            interval.file=gc.gene.file, output.file)
    }
    gatk.coverage <- output.file
}

### GC-normalize coverage -----------------------------------------------------

if (!is.null(gatk.coverage)) {
    library(PureCN)
    output.file <- file.path(outdir,  gsub(".txt$|_interval_summary",
        paste0("_", tolower(method), ".txt"), basename(gatk.coverage)))
    outpng.file <- sub("txt$","png", output.file)
    if (file.exists(output.file) && !force) {
        message(output.file, " exists. Skipping... (--force will overwrite)")
    } else {
        png(outpng.file, width=800)
        correctCoverageBias(gatk.coverage, gc.gene.file,
            output.file=output.file, method=method, plot.gc.bias=TRUE)
        dev.off()
   } 
}
    
