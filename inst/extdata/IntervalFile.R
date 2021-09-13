suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))

### Parsing command line ------------------------------------------------------

option_list <- list(
    make_option(c("--fasta"), action = "store", type = "character", 
        default = NULL, help = "Reference Fasta file"),
    make_option(c("--in-file"), action = "store", type = "character", 
        default = NULL,
        help = "Infile specifying target (baits) intervals. Needs to be parsable by rtracklayer."),
    make_option(c("--off-target"), action = "store_true",
        default = formals(PureCN::preprocessIntervals)$off.target, 
        help = "Include off-target regions [default %default]"),
    make_option(c("--average-target-width"), action = "store", type = "integer",
        default = formals(PureCN::preprocessIntervals)$average.target.width, 
        help = "Split large targets to approximately that size [default %default]"),
    make_option(c("--min-target-width"), action = "store", type = "integer",
        default = formals(PureCN::preprocessIntervals)$min.target.width, 
        help = "Either resize or drop targets smaller than specified [default %default]"),
    make_option(c("--small-targets"), action = "store", type = "character",
        default = formals(PureCN::preprocessIntervals)$small.targets[[2]], 
        help = "Either 'resize' or 'drop' small targets [default %default]"),
    make_option(c("--average-off-target-width"), action = "store", type = "integer",
        default = formals(PureCN::preprocessIntervals)$average.off.target.width, 
        help = "Bin off-target regions to approximately that size [default %default]"),
    make_option(c("--off-target-seqlevels"), action = "store", type = "character",
        default = formals(PureCN::preprocessIntervals)$off.target.seqlevels[[2]], 
        help = "Controls how to deal with chromosomes/contigs not found in --in-file. One of targeted, all [default %default]"),
    make_option(c("--mappability"), action = "store", type = "character", 
        help = "File parsable by rtracklayer specifying mappability scores of genomic regions."),
    make_option(c("--min-mappability"), action = "store", type = "character", 
        default = paste(eval(formals(PureCN::preprocessIntervals)$min.mappability), collapse=","),
        help = "Minimum mappability for on-target, off-target and chrY regions [default %default]"),
    make_option(c("--reptiming"), action = "store", type = "character", default = NULL,
        help = "File parsable by rtracklayer specifying replication timing scores of genomic regions."),
    make_option(c("--average-reptiming-width"), action = "store", type = "integer",
        default = formals(PureCN::preprocessIntervals)$average.reptiming.width, 
        help = "Average the replication timing data into bins of the specified size [default %default]"),
    make_option(c("--genome"), action = "store", type = "character", 
        default = NULL,
        help = "Genome version. If one of hg18, hg19, hg38, mm9, mm10, rn4, rn5, rn6, canFam3 will annotate intervals with gene symbols"),
    make_option(c("--out-file"), action = "store", type = "character", 
        default = NULL,
        help = "Outfile of annotated targets optimized for copy number calling."),
    make_option(c("--export"), action = "store", type = "character", default = NULL,
        help = "Export optimized intervals using rtracklayer. The file extension specifies the format."),
    make_option(c("--exclude"), action = "store", type = "character", 
        help = "File parsable by rtracklayer specifying baits that should be excluded from --in-file."),
    make_option(c("-v", "--version"), action = "store_true", default = FALSE,
        help = "Print PureCN version"),
    make_option(c("-f", "--force"), action = "store_true", default = FALSE,
        help = "Overwrite existing files")
)

alias_list <- list(
    "infile" = "in-file",
    "offtarget" = "off-target",
    "targetwidth" = "average-target-width",
    "mintargetwidth" = "min-target-width",
    "smalltargets" = "small-targets",
    "offtargetwidth" = "average-off-target-width",
    "offtargetseqlevels" = "off-target-seqlevels",
    "minmappability" = "min-mappability",
    "reptimingwidth" = "average-reptiming-width",
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
    q(status = 1)
}

outfile <- opt$out_file

if (is.null(opt$in_file)) stop("Need --in-file.")
if (is.null(opt$fasta)) stop("Need --fasta.")
if (is.null(outfile)) stop("Need --out-file.")

if (!opt$force && file.exists(outfile)) {
    stop(outfile, " exists. Use --force to overwrite.")
}

.checkOutputDir <-function(filename) {
    if (!is.null(filename)) {
        dname <- dirname(filename)
        if (!file.exists(dname)) {
            flog.fatal("Output directory %s does not exist.", dname)
            q(status = 1)
        }
    }
}

.checkOutputDir(outfile)
.checkOutputDir(opt$export)
.checkOutputDir(opt$exclude)

in.file <- normalizePath(opt$in_file, mustWork = TRUE)
reference.file <- normalizePath(opt$fasta, mustWork = TRUE)

suppressPackageStartupMessages(library(rtracklayer))

intervals <- try(import(in.file), silent = TRUE)

if (is(intervals, "try-error")) { 
    flog.warn("Could not parse --in-file with rtracklayer:\n\n%s\nTrying GATK3 parser that will probably fail...", intervals)
    intervals <- in.file
} else {
    seqlevels(intervals) <- seqlevelsInUse(intervals)
    if (sum(c("MT", "chrM", "chMT") %in% seqlevels(intervals))) {
        flog.warn("--in-file contains mitochondrion sequence. It is highly recommended to exclude those baits.")
    }
}    
mappability <- opt$mappability

if (!is.null(mappability)) {
    mappability <- normalizePath(mappability, mustWork = TRUE)
    flog.info("Loading %s...", mappability)
    mappability <- import(mappability)
}
exclude <- opt$exclude
if (!is.null(exclude)) {
    exclude <- normalizePath(exclude, mustWork = TRUE)
    flog.info("Loading %s...", exclude)
    exclude <- import(exclude)
}

flog.info("Loading PureCN %s...", Biobase::package.version("PureCN"))
suppressPackageStartupMessages(library(PureCN))
flog.info("Processing %s...", in.file)

seqinfoRef <- seqinfo(scanFaIndex(reference.file))

reptiming <- opt[["reptiming"]]
if (!is.null(reptiming)) {
    reptiming <- normalizePath(reptiming, mustWork = TRUE)
    flog.info("Loading %s...", reptiming)
    reptiming <- import(reptiming)
}

if (!opt$off_target) {
    flog.info("Will not add off-target regions. This is only recommended for%s",
     " Amplicon data. Add --off-target to include them.")
}

min.mappability <- as.numeric(strsplit(opt$min_mappability, ",")[[1]])

outGC <- preprocessIntervals(intervals, reference.file,
    output.file = outfile, off.target = opt$off_target,
    mappability = mappability, min.mappability = min.mappability,
    min.target.width = opt$min_target_width,
    small.targets = opt$small_targets,
    average.off.target.width = opt$average_off_target_width,
    reptiming = reptiming, average.reptiming.width = opt$average_reptiming_width,
    off.target.seqlevels = opt$off_target_seqlevels,
    exclude = exclude, average.target.width = opt$average_target_width)

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

.writeGc <- function(interval.gr, output.file) {
    tmp <- data.frame(
        Target = as.character(interval.gr),
        gc_bias = interval.gr$gc_bias,
        mappability = interval.gr$mappability,
        reptiming = interval.gr$reptiming,
        Gene = interval.gr$Gene,
        on_target = interval.gr$on.target
    )    
    write.table(tmp, file = output.file, row.names = FALSE, quote = FALSE, 
                sep = "\t")
}

if (!is.null(opt$genome) ) {
    if (is.null(knownGenome[[opt$genome]])) {
        flog.warn("%s genome not known. %s Known genomes: %s", opt$genome, 
        "Will not annotate targets with gene symbols.", 
        paste(names(knownGenome), collapse=", "))
    } else if (!require(knownGenome[[opt$genome]], character.only = TRUE)) {
        flog.warn("Install %s to get gene symbol annotation.", 
            knownGenome[[opt$genome]])
    } else if (!require(knownOrg[[opt$genome]], character.only = TRUE)) {
        flog.warn("Install %s to get gene symbol annotation.",
            knownOrg[[opt$genome]])
    } else {
        outGC <- suppressMessages(annotateTargets(outGC,
            get(knownGenome[[opt$genome]]), get(knownOrg[[opt$genome]])))
        .writeGc(outGC, outfile)
    }
} else {
    flog.warn("Specify --genome to get gene symbol annotation.")
}

if (!is.null(opt$export)) {
    export(outGC, opt$export)
}
