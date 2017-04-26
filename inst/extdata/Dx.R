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


# Parse input rds
infileRds <- opt$rds
if (is.null(infileRds)) stop("Need --rds")
infileRds <- normalizePath(infileRds, mustWork=TRUE)

# Parse outdir
outdir <- opt$outdir
if (is.null(outdir)) outdir <- dirname(infileRds)
outdir <- normalizePath(outdir, mustWork=TRUE)

# Parse both BED files restricting covered region
callableFile <- opt$callable
callable <- NULL
if (!is.null(callableFile)) {
    suppressPackageStartupMessages(library(rtracklayer))
    callableFile <- normalizePath(callableFile, mustWork=TRUE)
    flog.info("Reading %s...", callableFile)
    callable <- import(callableFile)
}

excludeFile <- opt$exclude
exclude <- NULL
if (!is.null(excludeFile)) {
    suppressPackageStartupMessages(library(rtracklayer))
    excludeFile <- normalizePath(excludeFile, mustWork=TRUE)
    flog.info("Reading %s...", excludeFile)
    exclude <- import(excludeFile)
}


### Run PureCN ----------------------------------------------------------------

flog.info("Loading PureCN...")
suppressPackageStartupMessages(library(PureCN))

flog.info("Reading %s...", infileRds)
res <- readCurationFile(infileRds)
sampleid <- res$input$sampleid

outfileMb <- file.path(outdir, paste0(sampleid, '_purecn_mutation_burden.csv'))

force <- !is.null(opt$force)

if (!force && file.exists(outfileMb)) {
    stop(outfileMb, " exists. Use --force to overwrite.")
}    

flog.info("Calling mutation burden...")
mb <- callMutationBurden(res, callable=callable, exclude=exclude)
write.csv(cbind(Sampleid=sampleid, mb), file=outfileMb, row.names=FALSE, quote=FALSE)
