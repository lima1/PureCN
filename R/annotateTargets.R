#' Annotate targets with gene symbols
#' 
#' This function can be used to add a \sQuote{Gene} meta column containing 
#' gene symbols to a \code{GRanges} object.
#' It applies heuristics to find the protein coding genes that were
#' likely meant to target in the assay design in case transcripts
#' overlap. 
#' 
#' @param x A \code{GRanges} object with interals to annotate
#' @param txdb A \code{TxDb} database, e.g. 
#' \code{TxDb.Hsapiens.UCSC.hg19.knownGene}
#' @param org A \code{OrgDb} object, e.g. \code{org.Hs.eg.db}.
#' @return A \code{GRanges} object.
#' @author Markus Riester
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(org.Hs.eg.db)
#' 
#' normal.coverage.file <- system.file("extdata", "example_normal.txt", 
#'     package="PureCN")
#' x <- head(readCoverageFile(normal.coverage.file),100)
#' x <- annotateTargets(x,TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db)
#' 
#' @importFrom GenomicFeatures transcriptsByOverlaps exonsByOverlaps
#' @export annotateTargets
annotateTargets <- function(x, txdb, org) {
    if (!is.null(x$on.target)) {
        idx <- x$on.target
    } else {
        idx <- seq_along(x)        
    }    
    txdb <- .checkSeqlevelStyle(x, txdb, "txdb", "interval file") 
    id <- transcriptsByOverlaps(txdb, ranges=x[idx], columns = "GENEID")
    id$SYMBOL <-suppressWarnings(select(org, sapply(id$GENEID, function(x)x[1]), "SYMBOL")[,2])
    idExons <- exonsByOverlaps(txdb, ranges=x[idx], columns = "GENEID")
    idExons$SYMBOL <-suppressWarnings(select(org, sapply(idExons$GENEID, function(x)x[1]), "SYMBOL")[,2])
	ov <- findOverlaps(x[idx], id)
	ovExons <- findOverlaps(x[idx], idExons)
    
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
    if (is.null(x$Gene)) x$Gene <- "."
    x[idx]$Gene[d.f$i] <- as.character(d.f$SYMBOL)
    x$Gene[is.na(x$Gene)] <- "."

    flog.warn("Attempted adding gene symbols to intervals. Heuristics have been %s",
        "used to pick symbols for overlapping genes.")
    x
}
