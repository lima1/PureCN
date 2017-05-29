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
'offtarget',    't', 0, "logical",
'accessible',   'b', 1, "character",
'genome',       'g', 1, "character"
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

intervals <- try(import(in.file), silent=TRUE)
if (class(intervals) == "try-error") intervals <- in.file

accessible <- opt$accessible

if (!is.null(accessible)) {
    accessible <- normalizePath(accessible, mustWork=TRUE)
    flog.info("Loading %s...", accessible)
    accessible <- import(accessible)
}
    
flog.info("Loading PureCN...")
suppressPackageStartupMessages(library(PureCN))
flog.info("Processing %s...", in.file)

if (is.null(opt$offtarget)) {
    flog.info("Will not add off-target regions. This is only recommended for%s",
     " Amplicon data. Add --offtarget to include them.")
}

outGC <- calculateGCContentByInterval(intervals, reference.file, 
    output.file = outfile, off.target=!is.null(opt$offtarget), 
    accessible=accessible)

knownGenome <- list(
    hg18="TxDb.Hsapiens.UCSC.hg18.knownGene",
    hg19="TxDb.Hsapiens.UCSC.hg19.knownGene",
    hg38="TxDb.Hsapiens.UCSC.hg38.knownGene",
    mm9="TxDb.Mmusculus.UCSC.mm9.knownGene",
    mm10="TxDb.Mmusculus.UCSC.mm10.knownGene"
)
knownOrg <- list(
    hg18="org.Hs.eg.db", 
    hg19="org.Hs.eg.db", 
    hg38="org.Hs.eg.db", 
    mm9="org.Mm.eg.db", 
    mm10="org.Mm.eg.db"
)

.writeGc <- function(interval.gr, output.file) {
    tmp <- data.frame(
        Target=as.character(interval.gr),
        gc_bias=interval.gr$gc_bias,
        Gene=interval.gr$Gene,
        on_target=interval.gr$on.target
    )    
    write.table(tmp, file=output.file, row.names=FALSE, quote=FALSE, sep="\t")
}

.annotateIntervals <- function(outGC, txdb, org, output.file = NULL) {
    idx <- outGC$on.target
    id <- transcriptsByOverlaps(txdb, ranges=outGC[idx], columns = "GENEID")
    id$SYMBOL <-suppressWarnings(select(org, sapply(id$GENEID, function(x)x[1]), "SYMBOL")[,2])
    outGC[idx]$Gene <- id$SYMBOL[findOverlaps(outGC[idx], id, select="first")]
    outGC$Gene[is.na(outGC$Gene)] <- "."
    .writeGc(outGC, output.file)
}
if (!is.null(opt$genome) ) {
    if (is.null(knownGenome[[opt$genome]])) {
        flog.warn("%s genome not known. %s", genome, 
        "Will not annotate targets with gene symbols.")
        q(status=1)
    }    
    if (!require(knownGenome[[opt$genome]], character.only=TRUE)) {
        flog.warn("Install %s to get gene symbol annotation.", 
            knownGenome[[opt$genome]])
        q(status=1)
    }
    if (!require(knownOrg[[opt$genome]], character.only=TRUE)) {
        flog.warn("Install %s to get gene symbol annotation.", 
            knownOrg[[opt$genome]])
        q(status=1)
    }
    .annotateIntervals(outGC, get(knownGenome[[opt$genome]]),
        get(knownOrg[[opt$genome]]), outfile)
} else {
    flog.warn("Specify --genome to get gene symbol annotation.")
}    
