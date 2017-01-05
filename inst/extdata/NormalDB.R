library('getopt')

### Parsing command line ------------------------------------------------------

spec <- matrix(c(
'help' ,  'h', 0, "logical",
'force' , 'f', 0, "logical",
'gcgene', 'c', 1, "character",
'method', 'm', 1, "character",
'coveragefiles', 'v', 1, "character",
'assay',   'a',1, "character",
'outdir' , 'o', 1, "character"
), byrow=TRUE, ncol=4)
opt <- getopt(spec)

if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE))
    q(status=1)
}

.checkFileList <- function(file) {
    files <- read.delim(file, as.is=TRUE, header=FALSE)[,1]
    numExists <- sum(file.exists(files), na.rm=TRUE)
    if (numExists < length(files)) { 
        stop("File not exists in file ", file)
    }
    files
}

if (is.null(opt$coveragefiles)) {
    stop("need --coveragefiles.")
}

coverageFiles <- .checkFileList(opt$coveragefiles)
outdir <- opt$outdir
if (is.null(outdir)) {
    stop("need --outdir")
}
outdir <- normalizePath(outdir, mustWork=TRUE)
method <- opt$method
assay <- opt$assay
if (is.null(assay)) assay <- ""

.getFileName <- function(outdir, prefix, suffix, assay, method) {
    if (nchar(assay)) assay <- paste0("_", assay)
    if (nchar(method)) method <- paste0("_", method)
    file.path(outdir, paste0(prefix, assay, method, suffix))
}

.gcNormalize <- function(gatk.coverage, gc.gene.file, method, outdir, force) {
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
     
if (length(coverageFiles)) {
    library(PureCN)
    if (!is.null(method)) {
        if (!method %in% c("LOESS", "POLYNOMIAL")) {
            stop("Unknown GC-normalization method")
        }
        gc.gene.file <- opt$gcgene
        gc.gene.file <- normalizePath(gc.gene.file, mustWork=TRUE)
        for (gatk.coverage in coverageFiles) {
            .gcNormalize(gatk.coverage, gc.gene.file, method, outdir, force)
        }    
    } else {
        method <- ""
        message("--method not specified, assuming coverage files are ",
            "GC-normalized")
    }        
    normalDB <- createNormalDatabase(coverageFiles)
    saveRDS(normalDB, file=.getFileName(outdir,"normalDB",".rds", assay, method))
}

if (length(coverageFiles) > 3) {
    library(PureCN)
    target.weight.file <- .getFileName(outdir,"target_weights",".txt", assay, method)
    createTargetWeights(coverageFiles[1:2], coverageFiles[-(1:2)], 
        target.weight.file)
} else {
    message("Not enough coverage files for creating target_weights.txt")
}

