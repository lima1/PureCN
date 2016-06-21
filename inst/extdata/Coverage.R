library('getopt')

spec <- matrix(c(
'help' , 'h', 0, "logical",
'bam', 'b', 1, "character",
'gatkcoverage', 'g', 1, "character",
'gcgene' , 'c', 1, "character",
'outdir' , 'o', 1, "character"
), byrow=TRUE, ncol=4)
opt <- getopt(spec)


if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE))
    q(status=1)
}

bam.file <- opt$bam
interval.file <- opt$interval
gatk.coverage <- opt$gatkcoverage
gc.gene.file <- opt$gcgene
outdir <- opt$outdir


outdir <- normalizePath(outdir, mustWork=TRUE)
gc.gene.file <- normalizePath(gc.gene.file, mustWork=TRUE)

if (!is.null(bam.file)) {
    library(PureCN)
    bam.file <- normalizePath(bam.file, mustWork=TRUE)
    output.file <- file.path(outdir,  gsub(".bam$","_coverage.txt", 
        basename(bam.file)))

    calculateBamCoverageByInterval(bam.file=bam.file, 
        interval.file=gc.gene.file, output.file)

    gatk.coverage <- output.file
}

if (!is.null(gatk.coverage)) {
    library(PureCN)
    output.file <- file.path(outdir,  gsub(".txt$","_loess.txt", 
        basename(gatk.coverage)))
    correctCoverageBias(gatk.coverage, gc.gene.file,
        output.file)
}
    
