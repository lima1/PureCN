library('getopt')
library(futile.logger)

### Parsing command line ------------------------------------------------------

spec <- matrix(c(
'help' ,        'h', 0, "logical",
'version',      'v', 0, "logical",
'force' ,       'f', 0, "logical",
'fasta',        'a', 1, "character",
'infile',       'i', 1, "character",
'outfile',      'o', 1, "character",
'subdivide',    's', 0, "logical"
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

force <- !is.null(opt$force)
outfile <- opt$outfile

if (!force && file.exists(outfile)) {
    stop(outfile, " exists. Use --force to overwrite.")
}    

if (is.null(opt$infile)) stop("Need --infile.")
if (is.null(opt$fasta)) stop("Need --fasta.")

in.file <- normalizePath(opt$infile, mustWork=TRUE)
reference.file <- normalizePath(opt$fasta, mustWork=TRUE)


suppressPackageStartupMessages(library(rtracklayer))

intervals <- import(in.file)

if (!is.null(opt$subdivide)) {
    flog.info("Subdividing large targets.")
    suppressPackageStartupMessages(library(exomeCopy))
    intervals <- subdivideGRanges(intervals, mean(width(intervals)))
} 
flog.info("Loading PureCN...")
suppressPackageStartupMessages(library(PureCN))
calculateGCContentByInterval(intervals, reference.file, 
    output.file = outfile)

