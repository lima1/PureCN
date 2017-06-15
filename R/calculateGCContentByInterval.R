#' Calculates GC content by interval
#' 
#' Uses \code{scanFa} from the Rsamtools package to retrieve GC content of
#' intervals in a reference FASTA file. Can optimize intervals for copy
#' number calling by tiling long intervals and by including off-target regions.
#' This optimization largely follows CNVkit. 
#' 
#' @param interval.file File specifying the intervals. Interval is expected in
#' first column in format CHR:START-END.  Instead of a file, a \code{GRanges}
#' object can be provided. This allows the use of BED files for example. Note
#' that GATK interval files are 1-based (first position of the genome is 1).
#' Other formats like BED files are often 0-based. The \code{import} function
#' will automatically convert to 1-based \code{GRanges}.
#' @param reference.file Reference FASTA file.
#' @param output.file Optionally, write GC content file.
#' @param off.target Include off-target regions.
#' @param average.target.width Split large targets to approximately this size.
#' @param min.off.target.width Only include off-target regions of that
#' size
#' @param average.off.target.width Split off-target regions to that
#' @param off.target.padding Pad off-target regions.
#' @param mappability Annotate intervals with mappability score. Assumed on a scale
#' from 0 to 1, with score being 1/(number alignments). Expected as \code{GRanges}
#' object with first meta column being the score. Regions outside these ranges are
#' ignored, assuming that \code{mappability} covers the whole accessible genome.
#' @param min.mappability \code{double(3)} specifying the minimum mappability score
#' for on-target, off-target, and chrY regions in that order. The chrY regions
#' are only used for sex determination in \sQuote{PureCN} and are therefore 
#' treated differently. Requires \code{mappability}.
#' @return Returns GC content by interval as \code{GRanges} object.
#' @author Markus Riester
#' @references Talevich et al. (2016). CNVkit: Genome-Wide Copy Number 
#' Detection and Visualization from Targeted DNA Sequencing. PLoS Comput Biol.
#'
#' @examples
#' 
#' reference.file <- system.file("extdata", "ex2_reference.fa", 
#'     package="PureCN", mustWork = TRUE)
#' interval.file <- system.file("extdata", "ex2_intervals.txt", 
#'     package="PureCN", mustWork = TRUE)
#' bed.file <- system.file("extdata", "ex2_intervals.bed", 
#'     package="PureCN", mustWork = TRUE)
#' calculateGCContentByInterval(interval.file, reference.file, 
#'     output.file="gc_file.txt")
#' 
#' intervals <- import(bed.file)
#' calculateGCContentByInterval(intervals, reference.file, 
#'     output.file="gc_file.txt")
#' 
#' @export calculateGCContentByInterval
#' @importFrom rtracklayer import
#' @importFrom Biostrings letterFrequency
#' @importFrom BiocGenerics unstrand
#' @importFrom stats aggregate
#' @importFrom S4Vectors mcols
calculateGCContentByInterval <- function(interval.file, reference.file,
output.file = NULL, off.target=FALSE, average.target.width=400, 
min.off.target.width=20000, average.off.target.width=200000,  
off.target.padding=-500, mappability=NULL, min.mappability=c(0.5,0.1,0.7)) {
    if (class(interval.file)=="GRanges") {
        interval.gr <- .checkIntervals(interval.file)
    } else {    
        interval.gr <- readCoverageFile(interval.file)
    }
    if (is.null(interval.gr$on.target)) interval.gr$on.target <- TRUE    
    
    containsOfftarget <- sum(interval.gr$on.target)!=length(interval.gr)

    if (containsOfftarget) {
        flog.warn("Intervals contain off-target regions. %s",
            "Will not change intervals.")
    } else {    
        # split large targets
        if (!is.null(average.target.width)) {
            tmp <- tile(interval.gr, width=average.target.width)
            interval.gr <- unlist(tmp)
            interval.gr$on.target <- TRUE
            nChanges <- sum(sapply(tmp, length)>1)
            if (nChanges > 0) {
                flog.info("Splitting %i large targets to an average width of %i.",
                    nChanges, average.target.width)
            }
        } 
        
        # find off-target regions
        if (off.target) {
            if (off.target.padding > 0) {
                .stopUserError("off.target.padding must be negative.")
            }    
            offRegions <- setdiff(scanFaIndex(reference.file), unstrand(interval.gr))
            if (!is.null(mappability)) {
                offRegions <- intersect(offRegions, mappability)
            }    
            offRegions <- offRegions[width(offRegions)>off.target.padding*-2]
            offRegions <- .padGranges(offRegions, off.target.padding)

            flog.info("Tiling off-target regions to an average width of %i.",
                average.off.target.width)
                
            offRegions <- unlist(tile(offRegions, width=average.off.target.width))
            offRegions$on.target <- FALSE
            offRegions <- offRegions[width(offRegions)>=min.off.target.width]
            offRegions <- offRegions[seqnames(offRegions) %in% seqlevels(interval.gr)]
            interval.gr <- merge(interval.gr, offRegions, all=TRUE, sort=TRUE)
        }    
    }

    interval.gr <- .annotateMappability(interval.gr, mappability, 
        min.mappability) 

    x <- scanFa(reference.file, interval.gr)
    GC.count <- letterFrequency(x,"GC")
    all.count <- letterFrequency(x,"ATGC")
    interval.gr$gc_bias <- as.vector(ifelse(all.count==0,NA,GC.count/all.count))
    # exclude unavailable regions
    interval.gr <- interval.gr[which(!is.na(interval.gr$gc_bias))]
    if (is.null(interval.gr$Gene)) interval.gr$Gene <- "."

    if (!is.null(output.file)) {
        .writeGc(interval.gr, output.file)
    }    
    invisible(interval.gr)
}


# add mappability score to intervals
.annotateMappability <- function(interval.gr, mappability, min.mappability) {    
    interval.gr$mappability <- 1
    if (!is.null(mappability)) {
        ov <- findOverlaps(interval.gr, mappability)
        mappScore <- aggregate(mcols(mappability)[subjectHits(ov),1], by=list(queryHits(ov)), mean)
        interval.gr$mappability[mappScore[,1]] <- mappScore[,2]
    } else {
        flog.warn("No mappability scores provided.")
        return(interval.gr)
    }    
    # remove intervals with low mappability
    nBefore <- sum(interval.gr$on.target)
    interval.gr <- interval.gr[
        (interval.gr$on.target & interval.gr$mappability >= min.mappability[1] ) |
        (!interval.gr$on.target & interval.gr$mappability >= min.mappability[2] ) ]
    # remove chrY low mappability    
    sex.chr <- .getSexChr(seqlevels(interval.gr))[2]
    interval.gr <- interval.gr[!seqnames(interval.gr) %in% sex.chr |
        interval.gr$mappability >= min.mappability[3] ]
    nAfter <- sum(interval.gr$on.target)
    if (nBefore > nAfter) {
        flog.info("Removing %i targets with low mappability score (<%.2f).", 
            nBefore-nAfter, min.mappability[1])
    }
    interval.gr
}

.writeGc <- function(interval.gr, output.file) {
    tmp <- data.frame(
        Target=as.character(interval.gr),
        gc_bias=interval.gr$gc_bias,
        mappability=interval.gr$mappability,
        Gene=interval.gr$Gene,
        on_target=interval.gr$on.target
    )    
    write.table(tmp, file=output.file, row.names=FALSE, quote=FALSE, sep="\t")
}
