library('getopt')
library(futile.logger)
library(BiocParallel)

### Parsing command line ------------------------------------------------------

spec <- matrix(c(
'help' ,        'h', 0, "logical",
'version',      'v', 0, "logical",
'force' ,       'f', 0, "logical",
'seed',         'S', 1, "integer", 
'cpu',          'C', 1, "integer", 
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
cpu <- if (is.null(opt$cpu)) 1 else opt$cpu

bam.file <- opt$bam
index.file <- opt$bai

gatk.coverage <- opt$gatkcoverage
gc.gene.file <- opt$gcgene
outdir <- opt$outdir
if (is.null(outdir)) outdir <- "."

method <- ifelse(is.null(opt$method), "LOESS", opt$method)

outdir <- normalizePath(outdir, mustWork=TRUE)
gc.gene.file <- normalizePath(gc.gene.file, mustWork=TRUE)

### Calculate coverage from BAM files -----------------------------------------

.checkFileList <- function(file) {
    files <- read.delim(file, as.is=TRUE, header=FALSE)[,1]
    numExists <- sum(file.exists(files), na.rm=TRUE)
    if (numExists < length(files)) { 
        stop("File not exists in file ", file)
    }
    files
}

getCoverageBams <- function(bamFiles, indexFiles, outdir, gc.gene.file, 
    force=FALSE, cpu=1) {

    bamFiles <- bamFiles
    indexFiles <- indexFiles
    outdir <- outdir
    gc.gene.file <- gc.gene.file
    force <- force

    if (cpu>1) flog.info("Using %i CPUs.", cpu)

    .getCoverageBam <- function(bam.file, index.file=NULL, outdir, 
        gc.gene.file, force) {
        output.file <- file.path(outdir,  gsub(".bam$","_coverage.txt", 
            basename(bam.file)))
        futile.logger::flog.info("Processing %s...", output.file)
        if (!is.null(index.file)) {
            index.file <- normalizePath(index.file, mustWork=TRUE)
            index.file <- sub(".bai$", "", index.file)
        } else if (file.exists(sub("bam$", "bai", bam.file))) {
            index.file <- sub(".bam$", "", bam.file)
        } else {    
            index.file <- bam.file
        }    
        #return(output.file)
        if (file.exists(output.file) && !force) {
            futile.logger::flog.info("%s exists. Skipping... (--force will overwrite)", output.file)
        } else {
            PureCN::calculateBamCoverageByInterval(bam.file=bam.file, 
                interval.file=gc.gene.file, output.file=output.file,
                index.file=index.file)
        }
        output.file
    }
    
    param <- new(class(bpparam()), workers=cpu)
    coverageFiles <- unlist(
        bplapply(seq_along(bamFiles), 
            function(i) .getCoverageBam(bamFiles[i], indexFiles[i], outdir, gc.gene.file, force), 
            BPPARAM=param)
    )

    coverageFiles
}

coverageFiles <- NULL
indexFiles <- NULL

flog.info("Loading PureCN...")
suppressPackageStartupMessages(library(PureCN))
    
if (!is.null(bam.file)) {
    bam.file <- normalizePath(bam.file, mustWork=TRUE)
    if (grepl(".list$", bam.file)) {
        bamFiles <- .checkFileList(bam.file)
        if (!is.null(index.file)) { 
            if(!grepl(".list$", index.file)) {
                stop("list of BAM files requires list of BAI files.")
            }
            indexFiles <- .checkFileList(index.file)   
        }
    } else {
        bamFiles <- bam.file
    }    
    if (length(bamFiles) != length(indexFiles) && !is.null(indexFiles)) {
        stop("List of BAM files and BAI files of different length.")
    }    

    coverageFiles <- getCoverageBams(bamFiles, indexFiles, outdir, 
        gc.gene.file, force, cpu) 
}

### GC-normalize coverage -----------------------------------------------------

.gcNormalize <- function(gatk.coverage, gc.gene.file, method, outdir, force) {
    output.file <- file.path(outdir,  gsub(".txt$|_interval_summary",
        paste0("_", tolower(method), ".txt"), basename(gatk.coverage)))
    outpng.file <- sub("txt$","png", output.file)
    if (file.exists(output.file) && !force) {
        flog.info("%s exists. Skipping... (--force will overwrite)", output.file)
    } else {
        png(outpng.file, width=800)
        correctCoverageBias(gatk.coverage, gc.gene.file,
            output.file=output.file, method=method, plot.gc.bias=TRUE)
        dev.off()
   } 
}

if (!is.null(gatk.coverage) || !is.null(coverageFiles)) {
    # started not from BAMs?
    if (is.null(coverageFiles)) {
        if (grepl(".list$", gatk.coverage)) {
            coverageFiles <- .checkFileList(gatk.coverage)
        } else {
            coverageFiles <- gatk.coverage
        }
    }
    for (gatk.coverage in coverageFiles)     
        .gcNormalize(gatk.coverage, gc.gene.file, method, outdir, force)
}
    
