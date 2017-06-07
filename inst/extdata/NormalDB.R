suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))

### Parsing command line ------------------------------------------------------

option_list <- list(
    make_option(c("-v", "--version"), action="store_true", default=FALSE, 
        help="Print PureCN version"),
    make_option(c("-f", "--force"), action="store_true", default=FALSE, 
        help="Overwrite existing files"),
    make_option(c("--maxmeancoverage"), action="store", type="integer", default=NULL,
        help="Maximum coverage (downscale samples exceeding this cutoff) [default auto]"),
    make_option(c("--coveragefiles"), action="store", type="character", default=NULL,
        help="List of input coverage files (supported formats: PureCN, GATK and CNVkit)"),
    make_option(c("--assay"), action="store", type="character", default="",
        help="Optional assay name used in output names [default %default]"),
    make_option(c("--genome"), action="store", type="character", default=NULL,
        help="Genome version, used in output names [default %default]"),
    make_option(c("--outdir"), action="store", type="character", default=NULL,
        help="Output directory to which results should be written")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (opt$version) {
    message(as.character(packageVersion("PureCN")))
    q(status=1)
}    

force <- !is.null(opt$force)

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
assay <- opt$assay
genome <- opt$genome
if (is.null(genome)) stop("Need --genome")

.getFileName <- function(outdir, prefix, suffix, assay, genome) {
    if (nchar(assay)) assay <- paste0("_", assay)
    if (nchar(genome)) genome <- paste0("_", genome)
    file.path(outdir, paste0(prefix, assay, genome, suffix))
}

flog.info("Loading PureCN...")
if (length(coverageFiles)) {
    suppressPackageStartupMessages(library(PureCN))
    flog.info("Creating normalDB. Assuming coverage files are GC-normalized.")
    normalDB <- createNormalDatabase(coverageFiles, max.mean.coverage=opt$maxmeancoverage)
    saveRDS(normalDB, file=.getFileName(outdir,"normalDB",".rds", assay, 
        genome))
}

if (length(coverageFiles) > 3) {
    suppressPackageStartupMessages(library(PureCN))
    target.weight.file <- .getFileName(outdir,"target_weights",".txt", assay, 
        genome)
    outpng.file <- sub("txt$","png", target.weight.file)
    flog.info("Creating target weights.")
    png(outpng.file, width=800, height=400)
    createTargetWeights(coverageFiles[1:2], coverageFiles[-(1:2)], 
        target.weight.file, plot=TRUE)
    dev.off()
} else {
    flog.warn("Not enough coverage files for creating target_weights.txt")
}
