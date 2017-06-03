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

if (is.null(opt$infile)) stop("Need --infile.")
if (is.null(opt$fasta)) stop("Need --fasta.")
if (is.null(opt$outfile)) stop("Need --outfile.")

if (!force && file.exists(outfile)) {
    stop(outfile, " exists. Use --force to overwrite.")
}    


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
    idExons <- exonsByOverlaps(txdb, ranges=outGC[idx], columns = "GENEID")
    idExons$SYMBOL <-suppressWarnings(select(org, sapply(idExons$GENEID, function(x)x[1]), "SYMBOL")[,2])
	ov <- findOverlaps(outGC[idx], id)
	ovExons <- findOverlaps(outGC[idx], idExons)
    
    # for targets with multiple gene hits, use the one with most overlapping
    # targets
	d.f <- data.frame(i=queryHits(ov), SYMBOL=as.character(id$SYMBOL[subjectHits(ov)]))
	d.f <- d.f[!duplicated(d.f),]

    # remove non-coding transcripts and ORFs
    d.f <- d.f[!grepl("-AS\\d$", d.f$SYMBOL),]
    d.f <- d.f[!grepl("^LOC\\d", d.f$SYMBOL),]
    d.f <- d.f[!grepl("\\dorf\\d", d.f$SYMBOL),]
    d.f <- d.f[!grepl("^FLJ\\d+$", d.f$SYMBOL),]

    d.f$COUNT <- table(d.f$SYMBOL)[d.f$SYMBOL]

    # in case multiple symbols have the same number of targets, prioritize the ones overlapping exons
    d.fExons <- data.frame(i = queryHits(ovExons), SYMBOL = as.character(idExons$SYMBOL[subjectHits(ovExons)]))
    d.f$OverlapsExon <- ifelse(paste(d.f[,1], d.f[,2]) %in% paste(d.fExons[,1], d.fExons[,2]), 1, 2)
    
    # reorder and pick the best transcript
	d.f <- d.f[order(d.f$i, d.f$COUNT, d.f$OverlapsExon),]
    d.f$FLAG <- duplicated(d.f$i, fromLast=TRUE)
	d.f <- d.f[!duplicated(d.f$i),]

    # Exclude targets for which we have multiple hits, but only one interval
    d.f <- d.f[!d.f$FLAG | d.f$COUNT>2,]
    outGC[idx]$Gene[d.f$i] <- as.character(d.f$SYMBOL)
    outGC$Gene[is.na(outGC$Gene)] <- "."

    flog.warn("Attempted adding gene symbols to intervals. Heuristics have been %s",
        "used to pick symbols for overlapping genes.")
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
