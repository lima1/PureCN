library('getopt')
library(futile.logger)

### Parsing command line ------------------------------------------------------

spec <- matrix(c(
'help' ,        'h', 0, "logical",
'version',      'v', 0, "logical",
'force' ,       'f', 0, "logical",
'rds',          'r', 1, "character",
'callable',     'a', 1, "character",
'exclude',      'b', 1, "character",
'out',          'o', 1, "character"
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


# Parse input rds
infileRds <- opt$rds
if (is.null(infileRds)) stop("Need --rds")
infileRds <- normalizePath(infileRds, mustWork=TRUE)

# Parse outdir

# Parse both BED files restricting covered region
callableFile <- opt$callable
callable <- NULL
if (!is.null(callableFile)) {
    suppressPackageStartupMessages(library(rtracklayer))
    callableFile <- normalizePath(callableFile, mustWork=TRUE)
    flog.info("Reading %s...", callableFile)
    callable <- GRanges(import(callableFile))
}

excludeFile <- opt$exclude
exclude <- NULL
if (!is.null(excludeFile)) {
    suppressPackageStartupMessages(library(rtracklayer))
    excludeFile <- normalizePath(excludeFile, mustWork=TRUE)
    flog.info("Reading %s...", excludeFile)
    exclude <- GRanges(import(excludeFile))
}


### Run PureCN ----------------------------------------------------------------

flog.info("Loading PureCN...")
suppressPackageStartupMessages(library(PureCN))

res <- readCurationFile(infileRds)
sampleid <- res$input$sampleid

.getOutPrefix <- function(opt, infile, sampleid) {
    out <- opt[["out"]]
    if (is.null(out)) {
        if (!is.null(infile) && file.exists(infile)) {
            outdir <- dirname(infile)
            prefix <- sampleid
        } else {
            stop("Need --out")
        }    
    } else {
        outdir <- dirname(out)
        prefix <- basename(out)
    }    
    outdir <- normalizePath(outdir, mustWork=TRUE)
    outPrefix <- file.path(outdir, prefix)
}
outPrefix <- .getOutPrefix(opt, infileRds, sampleid)
     
outfileMb <- paste0(outPrefix, '_mutation_burden.csv')

force <- !is.null(opt$force)

if (!force && file.exists(outfileMb)) {
    stop(outfileMb, " exists. Use --force to overwrite.")
}    

flog.info("Calling mutation burden...")
mb <- callMutationBurden(res, callable=callable, exclude=exclude)
write.csv(cbind(Sampleid=sampleid, mb), file=outfileMb, row.names=FALSE, quote=FALSE)
