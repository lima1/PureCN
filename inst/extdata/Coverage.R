library('getopt')
library(futile.logger)

### Parsing command line ------------------------------------------------------

spec <- matrix(c(
'help' ,        'h', 0, "logical",
'version',      'v', 0, "logical",
'force' ,       'f', 0, "logical",
'seed',         'S', 1, "integer", 
'bam',          'b', 1, "character",
'bai',          'a', 1, "character",
'gatkcoverage', 'g', 1, "character",
'gcgene',       'c', 1, "character",
'method',       'm', 1, "character",
'outdir',       'o', 1, "character"
), byrow=TRUE, ncol=4)
opt <- getopt(spec)

if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE))
    q(status=1)
}

if (!is.null(opt$version)) {
    message(as.character(packageVersion("PureCN")))
    q(status=1)
}    

if (!is.null(opt$seed)) {
    set.seed(opt$seed)
}
    

force <- !is.null(opt$force)

bam.file <- opt$bam
index.file <- opt$bai


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

    if (!is.null(index.file)) {
        index.file <- normalizePath(index.file, mustWork=TRUE)
        index.file <- sub(".bai$", "", index.file)
    } else if (file.exists(sub("bam$", "bai", bam.file))) {
        index.file <- sub(".bam$", "", bam.file)
    } else {    
        index.file <- bam.file
    }    

    if (file.exists(output.file) && !force) {
        flog.info("%s exists. Skipping... (--force will overwrite)", 
            output.file)
    } else {
        calculateBamCoverageByInterval(bam.file=bam.file, 
            interval.file=gc.gene.file, output.file=output.file,
             index.file=index.file)
    }
    gatk.coverage <- output.file
}

### GC-normalize coverage -----------------------------------------------------

.gcNormalize <- function(gatk.coverage, gc.gene.file, method, outdir, force) {
    output.file <- file.path(outdir,  gsub(".txt$|_interval_summary",
        paste0("_", tolower(method), ".txt"), basename(gatk.coverage)))
    outpng.file <- sub("txt$","png", output.file)
    if (file.exists(output.file) && !force) {
        flog.info("%s exists. Skipping... (--force will overwrite)", 
            output.file)
    } else {
        png(outpng.file, width=800)
        correctCoverageBias(gatk.coverage, gc.gene.file,
            output.file=output.file, method=method, plot.gc.bias=TRUE)
        dev.off()
   } 
}
.checkFileList <- function(file) {
    files <- read.delim(file, as.is=TRUE, header=FALSE)[,1]
    numExists <- sum(file.exists(files), na.rm=TRUE)
    if (numExists < length(files)) { 
        stop("File not exists in file ", file)
    }
    files
}

if (!is.null(gatk.coverage)) {
    library(PureCN)
    if (grepl(".list$", gatk.coverage)) {
        coverageFiles <- .checkFileList(gatk.coverage)
    } else {
        coverageFiles <- gatk.coverage
    }
    for (gatk.coverage in coverageFiles)     
        .gcNormalize(gatk.coverage, gc.gene.file, method, outdir, force)
}
    
