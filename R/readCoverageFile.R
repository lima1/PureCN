#' Read coverage file
#' 
#' Read coverage file produced by The Genome Analysis Toolkit or by
#' \code{\link{calculateBamCoverageByInterval}}.
#' 
#' 
#' @param file Target coverage file.
#' @param format File format. If missing, derived from the file 
#' extension. Currently only GATK DepthofCoverage format supported.
#' @return A \code{data.frame} with the parsed coverage information.
#' @author Markus Riester
#' @seealso \code{\link{calculateBamCoverageByInterval}}
#' @examples
#' 
#' tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
#'     package="PureCN")
#' coverage <- readCoverageFile(tumor.coverage.file)
#' 
#' @export readCoverageFile
readCoverageFile <- function(file, format) {
    if (missing(format)) format <- "GATK"
    targetCoverage <- .readCoverageGatk(file)
    .checkLowCoverage(targetCoverage)
    .checkIntervals(targetCoverage)
}


#' Read GATK coverage files
#' 
#' Read coverage file produced by The Genome Analysis Toolkit or by
#' \code{\link{calculateBamCoverageByInterval}}. This function
#' is deprecated. Please use \code{\link{readCoverageFile}} instead.
#' 
#' 
#' @param file Exon coverage file as produced by GATK.
#' @return A \code{data.frame} with the parsed coverage information.
#' @author Markus Riester
#' @seealso \code{\link{calculateBamCoverageByInterval}}
#' \code{\link{readCoverageFile}}
#' 
#' @export readCoverageGatk
readCoverageGatk <- function(file) {
    .Deprecated("readCoverageFile")
    readCoverageFile(file, format="GATK")
}

.readCoverageGatk <- function(file) {
    inputCoverage <- utils::read.table(file, header = TRUE)
    if (is.null(inputCoverage$total_coverage)) inputCoverage$total_coverage <- NA
    if (is.null(inputCoverage$average_coverage)) inputCoverage$average_coverage <- NA

    targetCoverage <- GRanges(inputCoverage$Target, 
        coverage=inputCoverage$total_coverage, 
        average.coverage=inputCoverage$average_coverage)

    targetCoverage
}

.checkIntervals <- function(coverageGr) {
    if (min(width(coverageGr))<2) {
        flog.warn("Coverage data contains single nucleotide intervals.")
    }    
    if(min(start(coverageGr))<1) {
        .stopUserError("Interval coordinates should start at 1, not at 0")
    }    
    ov <- findOverlaps(coverageGr, coverageGr)
    dups <- duplicated(queryHits(ov))
    
    if (sum(dups)) {
        dupLines <- subjectHits(ov)[dups]
        flog.warn("Found %i overlapping intervals, starting at line %i.",
            sum(dups), dupLines[1])
        coverageGr <- reduce(coverageGr)
    }
    coverageGr[order(coverageGr)]
}

.checkLowCoverage <- function(coverage) {
    chrsWithLowCoverage <- names(which(sapply(split(coverage$average.coverage, 
        as.character(seqnames(coverage))), mean, na.rm=TRUE) < 1))
    if (length(chrsWithLowCoverage)>2) {
        flog.warn("Multiple chromosomes with very low coverage: %s", 
            paste(chrsWithLowCoverage, collapse=","))
    }    
}

.addGCData <- function(tumor, gc.gene.file, verbose=TRUE) {
    tumor$gc_bias <- NA
    tumor$gc_bias <- as.numeric(tumor$gc_bias)
    tumor$Gene <- "."

    inputGC <- read.delim(gc.gene.file, as.is = TRUE)
    if (is.null(inputGC$gc_bias)) {
        .stopUserError("No gc_bias column in gc.gene.file.")
    }    
    if (is.null(inputGC$Gene)) {
        if (verbose) flog.info("No Gene column in gc.gene.file. You won't get gene-level calls.")
        inputGC$Gene <- "."
    }
    
    targetGC <- GRanges(inputGC[,1], gc_bias=inputGC$gc_bias, Gene=inputGC$Gene)

    ov <- findOverlaps(tumor, targetGC) 
    if (!identical(as.character(tumor), as.character(targetGC))) {
        # if only a few intervals are missing, maybe because of some poor 
        # quality regions, we just ignore those, otherwise we stop because 
        # user probably used the wrong file for the assay
        if (length(ov) < length(tumor)/2) {
            .stopUserError("tumor.coverage.file and gc.gene.file do not align.")
        } else {    
            flog.warn("tumor.coverage.file and gc.gene.file do not align.")
        }    
    }

    tumor[queryHits(ov)]$gc_bias <- targetGC[subjectHits(ov)]$gc_bias
    tumor[queryHits(ov)]$Gene <- targetGC[subjectHits(ov)]$Gene
    tumor <- .checkSymbolsChromosome(tumor)
    return(tumor)
}
    
