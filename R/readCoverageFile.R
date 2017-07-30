#' Read coverage file
#' 
#' Read coverage file produced by The Genome Analysis Toolkit or by
#' \code{\link{calculateBamCoverageByInterval}}.
#' 
#' 
#' @param file Target coverage file.
#' @param format File format. If missing, derived from the file 
#' extension. Currently only GATK DepthofCoverage format supported.
#' @param zero Start position is 0-based. Default is \code{FALSE}
#' for GATK, \code{TRUE} for BED file based intervals.
#' @return A \code{data.frame} with the parsed coverage information.
#' @author Markus Riester
#' @seealso \code{\link{calculateBamCoverageByInterval}}
#' @examples
#' 
#' tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
#'     package="PureCN")
#' coverage <- readCoverageFile(tumor.coverage.file)
#' 
#' @importFrom tools file_ext
#' @export readCoverageFile
readCoverageFile <- function(file, format, zero=NULL) {
    if (missing(format)) format <- .getFormat(file)
    if (format %in% c("cnn", "cnr")) {    
        targetCoverage <- .readCoverageCnn(file, zero, format)
    } else {
        targetCoverage <- .readCoverageGatk(file, zero)
    }
    .checkLowCoverage(targetCoverage)
    .checkIntervals(targetCoverage)
}

.getFormat <- function(file) {
    ext <- file_ext(file)
    if (ext %in% c("cnn", "cnr")) return(ext)
    "GATK"    
}    

#' Read GATK coverage files
#' 
#' This function is deprecated. Please use 
#' \code{\link{readCoverageFile}} instead.
#' 
#' @author Markus Riester
#' @seealso \code{\link{calculateBamCoverageByInterval}}
#' \code{\link{readCoverageFile}}
#' 
#' @export readCoverageGatk
readCoverageGatk <- function() {
    .Defunct("readCoverageFile")
}

.readCoverageGatk <- function(file, zero) {
    if (!is.null(zero)) flog.warn("zero ignored for GATK coverage files.")
    inputCoverage <- utils::read.table(file, header = TRUE)
    if (is.null(inputCoverage$total_coverage)) inputCoverage$total_coverage <- NA
    if (is.null(inputCoverage$counts)) inputCoverage$counts <- NA
    if (is.null(inputCoverage$on_target)) inputCoverage$on_target <- TRUE

    targetCoverage <- GRanges(inputCoverage$Target, 
        coverage=inputCoverage$total_coverage, 
        average.coverage=NA,
        counts=inputCoverage$counts,
        on.target=inputCoverage$on_target)
    targetCoverage <- .addAverageCoverage(targetCoverage)
    targetCoverage
}

.addAverageCoverage <- function(x) {
    x$average.coverage <- x$coverage/width(x)
    x
}    
.readCoverageCnn <- function(file, zero, format="cnn") {
    if (is.null(zero)) zero <- TRUE
    inputCoverage <- utils::read.table(file, header = TRUE)
    if (zero) inputCoverage$start <- inputCoverage$start + 1
    targetCoverage <- GRanges(inputCoverage)    
    targetCoverage$coverage <- targetCoverage$depth * width(targetCoverage)
    targetCoverage$average.coverage <- targetCoverage$depth
    targetCoverage$on.target <- TRUE
    targetCoverage$depth <- NULL
    targetCoverage$Gene <- targetCoverage$gene
    targetCoverage$on.target[which(targetCoverage$Gene=="Background")] <- FALSE
    targetCoverage$gene <- NULL
    if (format=="cnr") {
        targetCoverage$log.ratio <- targetCoverage$log2
    }
    targetCoverage$log2 <- NULL
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

    targets <- as.character(coverageGr)
    coverageGr <- sortSeqlevels(coverageGr)
    coverageGr <- sort(coverageGr)
    if (!identical(targets, as.character(coverageGr))) {
        flog.warn("Target intervals were not sorted.")
    }    
    # add fields that might be missing due to old PureCN versions
    if (is.null(coverageGr$counts)) coverageGr$counts <- NA
    if (is.null(coverageGr$on.target)) coverageGr$on.target <- TRUE
    coverageGr
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
    tumor$mappability <- 1
    tumor$gc_bias <- NA
    tumor$gc_bias <- as.numeric(tumor$gc_bias)
    if (is.null(tumor$Gene)) tumor$Gene <- "."

    inputGC <- read.delim(gc.gene.file, as.is = TRUE)
    if (is.null(inputGC$gc_bias)) {
        .stopUserError("No gc_bias column in gc.gene.file.")
    }    
    if (is.null(inputGC$Gene)) {
        if (verbose) flog.info("No Gene column in gc.gene.file. You won't get gene-level calls.")
        inputGC$Gene <- "."
    }
    if (is.null(inputGC$on_target)) {
        if (verbose) flog.info("No on_target column in gc.gene.file. Recreate this file with IntervalFile.R.")
        inputGC$on_target <- TRUE
    }
    if (is.null(inputGC$mappability)) {
        if (verbose) flog.info("No mappability column in gc.gene.file.")
        inputGC$mappability <- 1
    }
    
    targetGC <- GRanges(inputGC[,1], ranges=NULL, strand=NULL, inputGC[,-1])

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

    if(!is.null(tumor$on.target)) {
        if (!identical(tumor[queryHits(ov)]$on.target, targetGC[subjectHits(ov)]$on_target)) {
            flog.warn("Intervals in coverage and gc.gene.file have conflicting on/off-target annotation.")
            tumor[queryHits(ov)]$on.target <- targetGC[subjectHits(ov)]$on_target
        }
    } 
    tumor[queryHits(ov)]$mappability <- targetGC[subjectHits(ov)]$mappability
    tumor[queryHits(ov)]$gc_bias <- targetGC[subjectHits(ov)]$gc_bias
    tumor[queryHits(ov)]$Gene <- targetGC[subjectHits(ov)]$Gene
    tumor <- .checkSymbolsChromosome(tumor)
    return(tumor)
}
    
