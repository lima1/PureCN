suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))

### Parsing command line ------------------------------------------------------

option_list <- list(
    make_option(c("--fasta"), action = "store", type = "character", 
        default = NULL, help = "Reference Fasta file"),
    make_option(c("--infile"), action = "store", type = "character", 
        default = NULL,
        help = "Infile specifying target (baits) intervals. Needs to be parsable by rtracklayer."),
    make_option(c("--offtarget"), action = "store_true",
        default = formals(PureCN::preprocessIntervals)$off.target, 
        help = "Include off-target regions [default %default]"),
    make_option(c("--targetwidth"), action = "store", type = "integer",
        default = formals(PureCN::preprocessIntervals)$average.target.width, 
        help = "Split large targets to approximately that size [default %default]"),
    make_option(c("--offtargetwidth"), action = "store", type = "integer",
        default = formals(PureCN::preprocessIntervals)$average.off.target.width, 
        help = "Bin off-target regions to approximately that size [default %default]"),
    make_option(c("--offtargetseqlevels"), action = "store", type = "character",
        default = formals(PureCN::preprocessIntervals)$off.target.seqlevels[[2]], 
        help = "Controls how to deal with chromosomes/contigs not found in infile. One of targeted, all [default %default]"),
    make_option(c("--mappability"), action = "store", type = "character", 
        help = "File parsable by rtracklayer specifying mappability scores of genomic regions."),
    make_option(c("--minmappability"), action = "store", type = "character", 
        default = paste(eval(formals(PureCN::preprocessIntervals)$min.mappability), collapse=","),
        help = "Minimum mappability for on-target, off-target and chrY regions [default %default]"),
    make_option(c("--reptiming"), action = "store", type = "character", 
        help = "File parsable by rtracklayer specifying replication timing scores of genomic regions."),
    make_option(c("--reptimingbinsize"), action = "store", type = "integer", default = 100000,
        help = "Average the replication timing data into bins of the specified size [default %default]"),
    make_option(c("--genome"), action = "store", type = "character", 
        default = NULL,
        help = "Genome version. If one of hg18, hg19, hg38, mm9, mm10, rn4, rn5, rn6 will annotate intervals with gene symbols"),
    make_option(c("--outfile"), action = "store", type = "character", 
        default = NULL,
        help = "Outfile of annotated targets optimized for copy number calling."),
    make_option(c("--export"), action = "store", type = "character", default = NULL,
        help = "Export optimized intervals using rtracklayer. The file extension specifies the format."),
    make_option(c("--exclude"), action = "store", type = "character", 
        help = "File parsable by rtracklayer specifying baits that should be excluded from --infile."),
    make_option(c("-v", "--version"), action = "store_true", default = FALSE,
        help = "Print PureCN version"),
    make_option(c("-f", "--force"), action = "store_true", default = FALSE,
        help = "Overwrite existing files")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (opt$version) {
    message(as.character(packageVersion("PureCN")))
    q(status = 1)
}

outfile <- opt$outfile

if (is.null(opt$infile)) stop("Need --infile.")
if (is.null(opt$fasta)) stop("Need --fasta.")
if (is.null(opt$outfile)) stop("Need --outfile.")

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

.checkOutputDir(opt$outfile)
.checkOutputDir(opt$export)
.checkOutputDir(opt$exclude)

in.file <- normalizePath(opt$infile, mustWork = TRUE)
reference.file <- normalizePath(opt$fasta, mustWork = TRUE)

suppressPackageStartupMessages(library(rtracklayer))

intervals <- try(import(in.file), silent = TRUE)
if (class(intervals) == "try-error") { 
    intervals <- in.file
} else {
    if (sum(c("MT", "chrM", "chMT", "MT") %in% seqlevels(intervals))) {
        flog.warn("--infile contains mitochondrion sequence. It is highly recommended to exclude those baits.")
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
    if (opt$reptimingbinsize > 0) {
        flog.info("Averaging reptiming into bins of size %i...", opt$reptimingbinsize)
        bins <- tileGenome(seqinfoRef, tilewidth = opt$reptimingbinsize, 
            cut.last.tile.in.chrom=TRUE)
        reptiming <- PureCN:::.addScoreToGr(bins, reptiming, "score")
        reptiming <- reptiming[!is.na(score(reptiming))]
    }
}

if (!opt$offtarget) {
    flog.info("Will not add off-target regions. This is only recommended for%s",
     " Amplicon data. Add --offtarget to include them.")
}

min.mappability <- as.numeric(strsplit(opt$minmappability, ",")[[1]])

outGC <- preprocessIntervals(intervals, reference.file,
    output.file = outfile, off.target = opt$offtarget,
    mappability = mappability, min.mappability = min.mappability,
    average.off.target.width = opt$offtargetwidth,
    reptiming = reptiming, off.target.seqlevels = opt$offtargetseqlevels,
    exclude = exclude, average.target.width = opt$targetwidth)

knownGenome <- list(
    hg18 = "TxDb.Hsapiens.UCSC.hg18.knownGene",
    hg19 = "TxDb.Hsapiens.UCSC.hg19.knownGene",
    hg38 = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    mm9 = "TxDb.Mmusculus.UCSC.mm9.knownGene",
    mm10 = "TxDb.Mmusculus.UCSC.mm10.knownGene",
    rn4 = "TxDb.Rnorvegicus.UCSC.rn4.ensGene",
    rn5 = "TxDb.Rnorvegicus.UCSC.rn5.ensGene",
    rn6 = "TxDb.Rnorvegicus.UCSC.rn6.ensGene"
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
