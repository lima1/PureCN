#' Read allelic counts file
#' 
#' Read file containing counts of ref and alt alleles of common
#  SNPs by external tools like The Genome Analysis 
#' Toolkit 4. 
#' 
#' @param file Input file containing counts of ref and alt alleles
#' @param format File format. If missing, derived from the file 
#' extension. Currently only GATK4 CollectAllelicCounts (tsv)
#' format supported.
#' @param zero Start position is 0-based. Default is \code{FALSE}
#' for GATK, \code{TRUE} for BED file based intervals.
#' @return A \code{CollapsedVCF} with the parsed allelic counts.
#' @author Markus Riester
#' @examples
#' 
#' ac.file <- system.file("extdata", "example_allelic_counts.tsv", 
#'     package="PureCN")
#' vcf_ac <- readAllelicCountsFile(ac.file)
#' 
#' @importFrom utils write.table
#' @importFrom Biostrings DNAStringSet DNAStringSetList
#' @export readAllelicCountsFile
readAllelicCountsFile <- function(file, format, zero=NULL) {
    if (missing(format)) format <- "tsv"
    .readAllelicCountsFileGatk4(file, zero)
}

.writeAllelicCountsFileGatk <- function(vcf, id = 1, file) {
    outputCounts <- data.frame(
        CONTIG = seqnames(vcf),
        POSITION = start(vcf),
        REF_COUNT = sapply(geno(vcf)$AD[,id], function(x) x[1]),
        ALT_COUNT = sapply(geno(vcf)$AD[,id], function(x) x[2]),
        REF_NUCLEOTIDE = as.character(ref(vcf)),
        ALT_NUCLEOTIDE = unlist(CharacterList(alt(vcf)))
    )
    con <- file(file, open = "w")
    writeLines(paste("@RG", "ID:PureCN", paste0("SM:", samples(header(vcf))[id]), sep="\t"), con)
    write.table(outputCounts, con, row.names = FALSE, quote = FALSE, sep = "\t")
    close(con)
    invisible(outputCounts)
}

.readAllelicCountsFileGatk4 <- function(file, zero) {
    .extractSampleId <- function(line) {
        fields <- strsplit(line, "\t")[[1]]
        fields <- fields[grep("^SM:", fields)]
        gsub("^SM:", "", fields[1]) 
    }

    if (!is.null(zero)) flog.warn("zero ignored for GATK4 files.")
    con <- file(file, open = "r")
    sid <- NULL
    while ( TRUE ) {
        line = readLines(con, n = 1)
        if ( length(line) == 0 || !grepl("^@", line)[1]) {
            break
        }
        if (grepl("^@RG", line)[1]) sid <- .extractSampleId(line)
    }

    inputCounts <- read.delim(con, header = FALSE)
    colnames(inputCounts) <- strsplit(line, "\t")[[1]]
    close(con)
    inputCounts <-inputCounts[which(inputCounts$ALT_COUNT > 0),]
    gr <- GRanges(seqnames = inputCounts$CONTIG, IRanges(start = inputCounts$POSITION, end = inputCounts$POSITION))
    vcf <- VCF(gr, 
                colData = DataFrame(Samples = 1, row.names = sid),
                exptData = list(header = VCFHeader(samples = sid)))
    ref(vcf) <- DNAStringSet(inputCounts$REF_NUCLEOTIDE)
    alt(vcf) <- DNAStringSetList(lapply(inputCounts$ALT_NUCLEOTIDE, DNAStringSet))
    info(header(vcf)) <- DataFrame(
        Number ="0",
        Type = "Flag",
        Description = "dbSNP Membership",
        row.names="DB")

    geno(header(vcf)) <- DataFrame(
        Number =".",
        Type = "Integer",
        Description = "Allelic depths for the ref and alt alleles in the order listed",
        row.names="AD")

    info(vcf)$DB <- TRUE
    geno(vcf)$AD <- matrix(lapply(seq(nrow(inputCounts)), function(i) c( inputCounts$REF_COUNT[i], inputCounts$ALT_COUNT[i])), ncol = 1, dimnames = list(NULL, sid))
    names(vcf) <- paste0(seqnames(vcf), ":", start(vcf))
    .readAndCheckVcf(vcf)
}

