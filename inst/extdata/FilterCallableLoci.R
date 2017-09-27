# Simple helper script that only keeps CDS regions from an input BED file

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))

### Parsing command line ------------------------------------------------------

option_list <- list(
    make_option(c("--infile"), action="store", type="character", default=NULL,
        help="Infile specifying callable intervals. Needs to be parsable by rtracklayer."),
    make_option(c("--genome"), action="store", type="character", default=NULL,
        help="Genome version to filter non-CDS region. One of hg18, hg19, hg38, mm9, mm10, rn4, rn5, rn6."),
    make_option(c("--outfile"), action="store", type="character", default=NULL,
        help="Outfile with overlapping CDS regions in BED format."),
    make_option(c("-v", "--version"), action="store_true", default=FALSE, 
        help="Print PureCN version"),
    make_option(c("-f", "--force"), action="store_true", default=FALSE, 
        help="Overwrite existing files")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (opt$version) {
    message(as.character(packageVersion("PureCN")))
    q(status=1)
}    

outfile <- opt$outfile

if (is.null(opt$infile)) stop("Need --infile.")
if (is.null(opt$genome)) stop("Need --genome.")
if (is.null(opt$outfile)) stop("Need --outfile.")

if (!opt$force && file.exists(outfile)) {
    stop(outfile, " exists. Use --force to overwrite.")
}    

in.file <- normalizePath(opt$infile, mustWork=TRUE)

suppressPackageStartupMessages(library(rtracklayer))

intervals <- try(import(in.file), silent=TRUE)
if (class(intervals) == "try-error") intervals <- in.file

knownGenome <- list(
    hg18="TxDb.Hsapiens.UCSC.hg18.knownGene",
    hg19="TxDb.Hsapiens.UCSC.hg19.knownGene",
    hg38="TxDb.Hsapiens.UCSC.hg38.knownGene",
    mm9="TxDb.Mmusculus.UCSC.mm9.knownGene",
    mm10="TxDb.Mmusculus.UCSC.mm10.knownGene",
    rn4="TxDb.Rnorvegicus.UCSC.rn4.ensGene",
    rn5="TxDb.Rnorvegicus.UCSC.rn5.ensGene",
    rn6="TxDb.Rnorvegicus.UCSC.rn6.ensGene"
)

if (is.null(knownGenome[[opt$genome]])) {
    flog.warn("%s genome not known. %s", genome)
} else if (!require(knownGenome[[opt$genome]], character.only=TRUE)) {
    flog.warn("Install %s.", knownGenome[[opt$genome]])
} else {
    coding <- cds(get(knownGenome[[opt$genome]]))
    seqlevelsStyle(coding) <- seqlevelsStyle(intervals)
    export(reduce(intersect(intervals, unstrand(coding))), opt$outfile)
}

