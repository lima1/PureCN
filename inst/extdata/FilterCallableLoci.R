# Simple helper script that only keeps CDS regions from an input BED file

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))

### Parsing command line ------------------------------------------------------

option_list <- list(
    make_option(c("--in-file"), action = "store", type = "character", default = NULL,
        help = "Infile specifying callable intervals. Needs to be parsable by rtracklayer."),
    make_option(c("--genome"), action = "store", type = "character", default = NULL,
        help="Genome version to filter non-CDS region. One of hg18, hg19, hg38, mm9, mm10, rn4, rn5, rn6, canFam3."),
    make_option(c("--exclude"), action="store", type="character", default=NULL,
        help="Regular expression matching gene symbols to ignore, e.g. '^HLA'. Can also be a file parsable by rtracklayer."),
    make_option(c("--out-file"), action="store", type="character", default=NULL,
        help="Outfile with overlapping CDS regions in BED format."),
    make_option(c("-v", "--version"), action="store_true", default=FALSE, 
        help="Print PureCN version"),
    make_option(c("-f", "--force"), action="store_true", default=FALSE, 
        help="Overwrite existing files")
)

alias_list <- list(
    "infile" = "in-file",
    "outfile" = "out-file"
)
replace_alias <- function(x, deprecated = TRUE) {
    idx <- match(x, paste0("--", names(alias_list)))
    if (any(!is.na(idx))) {
        replaced <- paste0("--", alias_list[na.omit(idx)])
        x[!is.na(idx)] <- replaced
        if (deprecated) {
            flog.warn("Deprecated arguments, use %s instead.", paste(replaced, collapse=" "))
        }
    }
    return(x)
}
    
opt <- parse_args(OptionParser(option_list = option_list),
    args = replace_alias(commandArgs(trailingOnly = TRUE)),
    convert_hyphens_to_underscores = TRUE)

if (opt$version) {
    message(as.character(packageVersion("PureCN")))
    q(status = 0)
}


if (is.null(opt$in_file)) stop("Need --in-file.")
if (is.null(opt$genome)) stop("Need --genome.")
if (is.null(opt$out_file)) stop("Need --out-file.")

if (!opt$force && file.exists(opt$out_file)) {
    stop(opt$out_file, " exists. Use --force to overwrite.")
}    

in.file <- normalizePath(opt$in_file, mustWork = TRUE)

suppressPackageStartupMessages(library(rtracklayer))

intervals <- try(import(in.file), silent = TRUE)
if (is(intervals, "try-error")) {
    flog.warn("Could not parse --in-file with rtracklayer:\n\n%s\nTrying GATK3 parser that will probably fail...", intervals)
    intervals <- in.file
}
knownGenome <- list(
    hg18 = "TxDb.Hsapiens.UCSC.hg18.knownGene",
    hg19 = "TxDb.Hsapiens.UCSC.hg19.knownGene",
    hg38 = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    mm9 = "TxDb.Mmusculus.UCSC.mm9.knownGene",
    mm10 = "TxDb.Mmusculus.UCSC.mm10.knownGene",
    rn4 = "TxDb.Rnorvegicus.UCSC.rn4.ensGene",
    rn5 = "TxDb.Rnorvegicus.UCSC.rn5.ensGene",
    rn6 = "TxDb.Rnorvegicus.UCSC.rn6.ensGene",
    canFam3 = "TxDb.Cfamiliaris.UCSC.canFam3.refGene"
)

knownOrg <- list(
    hg18 = "org.Hs.eg.db",
    hg19 = "org.Hs.eg.db",
    hg38 = "org.Hs.eg.db", 
    mm9 = "org.Mm.eg.db",
    mm10 = "org.Mm.eg.db",
    rn4 = "org.Rn.eg.db",
    rn5 = "org.Rn.eg.db",
    rn6 = "org.Rn.eg.db",
    canFam3 = "org.Cf.eg.db"
)

flog.info("Loading %s...", knownGenome[[opt$genome]])
if (is.null(knownGenome[[opt$genome]])) {
    flog.warn("%s genome not known. %s", genome)
} else if (!suppressPackageStartupMessages(require(knownGenome[[opt$genome]], character.only = TRUE))) {
    flog.warn("Package %s not found.", knownGenome[[opt$genome]])
} else {
    coding <- cds(get(knownGenome[[opt$genome]]))
    seqlevelsStyle(coding) <- PureCN:::.getSeqlevelsStyle(intervals)
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
    export(intervalsCDS, opt$out_file)
    flog.info("Total size of CDS region: %.2fMb (%.2fMb input).", 
        PureCN:::.calcTargetedGenome(intervalsCDS), 
        PureCN:::.calcTargetedGenome(intervals))
}
