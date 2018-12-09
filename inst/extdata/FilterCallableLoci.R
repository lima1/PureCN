# Simple helper script that only keeps CDS regions from an input BED file

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))

### Parsing command line ------------------------------------------------------

option_list <- list(
    make_option(c("--infile"), action="store", type="character", default=NULL,
        help="Infile specifying callable intervals. Needs to be parsable by rtracklayer."),
    make_option(c("--genome"), action="store", type="character", default=NULL,
        help="Genome version to filter non-CDS region. One of hg18, hg19, hg38, mm9, mm10, rn4, rn5, rn6."),
    make_option(c("--exclude"), action="store", type="character", default=NULL,
        help="Regular expression matching gene symbols to ignore, e.g. '^HLA'. Can also be a file parsable by rtracklayer."),
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
if (is(intervals, "try-error")) intervals <- in.file

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
knownOrg <- list(
    hg18 = "org.Hs.eg.db",
    hg19 = "org.Hs.eg.db",
    hg38 = "org.Hs.eg.db", 
    mm9 = "org.Mm.eg.db",
    mm10 = "org.Mm.eg.db",
    rn4 = "org.Rn.eg.db",
    rn5 = "org.Rn.eg.db",
    rn6 = "org.Rn.eg.db"
)

flog.info("Loading %s...", knownGenome[[opt$genome]])
if (is.null(knownGenome[[opt$genome]])) {
    flog.warn("%s genome not known. %s", genome)
} else if (!suppressPackageStartupMessages(require(knownGenome[[opt$genome]], character.only=TRUE))) {
    flog.warn("Package %s not found.", knownGenome[[opt$genome]])
} else {
    coding <- cds(get(knownGenome[[opt$genome]]))
    seqlevelsStyle(coding) <- seqlevelsStyle(intervals)[1]
    intervalsCDS <- reduce(intersect(intervals, unstrand(coding)))

    if (!is.null(opt$exclude)) {
        if (!require(knownOrg[[opt$genome]], character.only = TRUE)) {
            flog.warn("Install %s to get gene symbol filtering.",
            knownOrg[[opt$genome]])
        } else {
            if (file.exists(opt$exclude)) {
                exclude <- import(opt$exclude)
            } else {     
                suppressPackageStartupMessages(library(PureCN))
                exclude <- suppressMessages(annotateTargets(intervalsCDS,
                    get(knownGenome[[opt$genome]]), get(knownOrg[[opt$genome]])))
                exclude <- reduce(exclude[grep(opt$exclude, exclude$Gene)])
            }
            flog.info("Excluded region of size %ibp.", sum(width(exclude)))
            intervalsCDS <- setdiff(intervalsCDS, exclude)
        }
    }    
    export(intervalsCDS, opt$outfile)
    flog.info("Total size of CDS region: %.2fMb (%.2fMb input).", 
        PureCN:::.calcTargetedGenome(intervalsCDS), 
        PureCN:::.calcTargetedGenome(intervals))
}

