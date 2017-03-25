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
    coverageData <- .readCoverageGatk(file)
    .checkLowCoverage(coverageData)
    return(.checkIntervals(coverageData))
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
    # originally from the ExomeCNV package
    gatk <- utils::read.table(file, header = TRUE)
    chrpos <- matrix(unlist(strsplit(as.character(gatk$Target), 
        ":")), ncol = 2, byrow = TRUE)
    
    idx <- !grepl("-", chrpos[,2])

    if (sum(idx)) {
        warning("Coverage data contains single nucleotide intervals.")
    }    
    chrpos[idx,2] <- paste0(chrpos[idx,2], "-", chrpos[idx,2])
    chr <- factor(chrpos[, 1])
    
    pos <- matrix(as.integer(unlist(strsplit(chrpos[, 2], "-"))), 
        ncol = 2, byrow = TRUE)
    start <- pos[, 1]
    end <- pos[, 2]

    if (is.null(gatk$total_coverage)) gatk$total_coverage <- NA
    if (is.null(gatk$average_coverage)) gatk$average_coverage <- NA

    tmp <- data.frame(probe = gatk$Target, chr = chr, probe_start = start, 
        probe_end = end, targeted.base = end - start + 1, sequenced.base = NA, 
        coverage = as.numeric(gatk$total_coverage), 
        average.coverage = as.numeric(gatk$average_coverage), 
        base.with..10.coverage = NA)
    return(tmp)
}


.checkIntervals <- function(coverage) {
    coverageGr <- GRanges(coverage$chr, IRanges(start=coverage$probe_start,
        end=coverage$probe_end))
    ov <- findOverlaps(coverageGr, coverageGr)
    dups <- duplicated(queryHits(ov))
    
    if (sum(dups)) {
        dupLines <- subjectHits(ov)[dups]
        wOv <- width(pintersect(coverageGr[queryHits(ov)], coverageGr[subjectHits(ov)]))
        singleDups <- which(wOv[dups]==1)
        sumSingle <- sum(wOv[dups]==1)
        tmpMsg <- ""
        if (sumSingle > 0) {
            tmpMsg <- paste(sumSingle, "intervals overlap by just 1bp.")
        }    
        warning("Found ", sum(dups), " overlapping intervals, starting ",
            "at line ", dupLines[1], ".\nKeeping them",
            " for now, but you should fix the interval file. ", 
            "Future PureCN versions might stop here. ",
            tmpMsg)
        #coverage <- coverage[-dupLines,]
    }
    coverage     
}

.checkLowCoverage <- function(coverage) {
    chrsWithLowCoverage <- names(which(sapply(split(coverage$average.coverage, 
        coverage$chr), mean, na.rm=TRUE) < 1))
    if (length(chrsWithLowCoverage)>2) {
        warning("Multiple chromosomes with very low coverage: ", 
            paste(chrsWithLowCoverage, collapse=","))
    }    
}
