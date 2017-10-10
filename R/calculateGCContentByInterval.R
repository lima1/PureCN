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
#' @param off.target.seqlevels Controls how to deal with chromosomes/contigs
#' found in the \code{reference.file} but not in the \code{interval.file}.
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
#' @importFrom GenomeInfoDb seqlevelsInUse seqlengths seqlevels<-
calculateGCContentByInterval <- function(interval.file, reference.file,
output.file = NULL, off.target=FALSE, average.target.width=400, 
min.off.target.width=20000, average.off.target.width=200000,  
off.target.padding=-500, mappability=NULL, min.mappability=c(0.5,0.1,0.7),
off.target.seqlevels=c("targeted", "noncircular", "all")) {
    if (class(interval.file)=="GRanges") {
        interval.gr <- .checkIntervals(interval.file)
    } else {    
        interval.gr <- readCoverageFile(interval.file)
    }
    if (is.null(interval.gr$on.target)) interval.gr$on.target <- TRUE    
    
    # make sure the chromsome naming style is the same in all 3 files
    # be nice and fix it if necessary
    interval.gr <- .checkSeqlevelStyle(scanFaIndex(reference.file), interval.gr, "interval")
    if (!is.null(mappability)) {
        mappability <- .checkSeqlevelStyle(scanFaIndex(reference.file), mappability, "mappability")
        mappability <- .remove0MappabilityRegions(mappability)
    }
     
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
            nChanges <- sum(elementNROWS(tmp) > 1)

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
            offRegions <- .dropShortUntargeted(offRegions, interval.gr)

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
            off.target.seqlevels <- match.arg(off.target.seqlevels)
            seqlevelsBefore <- seqlevelsInUse(offRegions)
            if (off.target.seqlevels == "targeted") {
                offRegions <- offRegions[seqnames(offRegions) %in% seqlevels(interval.gr)]
            } else if (off.target.seqlevels == "noncircular") {
                offRegions <- offRegions[seqnames(offRegions) %in% 
                    .getNonCircularSeqnames(reference.file)]
            }    
            seqlevelsAfter <- seqlevelsInUse(offRegions)
            if (!identical(seqlevelsBefore, seqlevelsAfter)) {
                flog.info("Removing following contigs from off-target regions: %s", 
                    paste(setdiff(seqlevelsBefore, seqlevelsAfter), collapse=","))
            }    
            interval.gr <- merge(offRegions, interval.gr, all = TRUE, sort=TRUE)
        }    
    }

    interval.gr <- .annotateMappability(interval.gr, mappability, 
        min.mappability) 
    
    flog.info("Calculating GC-content...")
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

# this function removes short chromosomes that have no probes (mainly a
# general way to remove chrM)
.dropShortUntargeted <- function(offRegions, interval.gr) {
    idx <- seqlevels(offRegions) %in% seqlevels(interval.gr) |
        seqlengths(offRegions) > 100000
    offRegions <- offRegions[seqnames(offRegions) %in% seqlevels(offRegions)[idx]]
    seqlevels(offRegions) <- seqlevelsInUse(offRegions)
    offRegions
}
    
# add mappability score to intervals
.annotateMappability <- function(interval.gr, mappability, min.mappability) {    
    interval.gr$mappability <- 1
    if (!is.null(mappability)) {
        ov <- findOverlaps(interval.gr, mappability)
        colScore <- if (is.null(mappability$score)) 1 else "score"
        mappScore <- aggregate(mcols(mappability)[subjectHits(ov),colScore], by=list(queryHits(ov)), mean)
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

.remove0MappabilityRegions <- function(mappability) {
    colScore <- if (is.null(mappability$score)) 1 else "score"
    mappability[which(mcols(mappability)[, colScore]>0),]
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

.checkSeqlevelStyle <- function(ref, x, name1, name2="reference") {
    refSeqlevelStyle <- try(seqlevelsStyle(ref), silent=TRUE)
    # if unknown, we cannot check and correct
    if (class(refSeqlevelStyle) == "try-error") return(x)
    xSeqlevelStyle <- try(seqlevelsStyle(x), silent=TRUE)

    if (class(xSeqlevelStyle) == "try-error") {
        .stopUserError("Chromosome naming style of ", name1, 
            " file unknown, should be ", refSeqlevelStyle, ".") 
    }    

    if (!length(intersect(xSeqlevelStyle,refSeqlevelStyle))) {
        flog.warn("Chromosome naming style of %s file (%s) was different from %s (%s).", 
            name1, xSeqlevelStyle, name2, refSeqlevelStyle)
        seqlevelsStyle(x) <- refSeqlevelStyle[1]
    }        
    x
}       

.getNonCircularSeqnames <- function(reference.file) {
    gs <- genomeStyles()
    style <- try(seqlevelsStyle(scanFaIndex(reference.file)), silent=TRUE)
    # style not known? then don't exclude any chromosomes, return all
    if (class(style) == "try-error") {
        return(seqlevels(scanFaIndex(reference.file)))
    }    
    gs <- gs[sapply(gs, function(x) style %in% colnames(x))]
    gss <- lapply(gs, function(x) list(seqlevels(scanFaIndex(reference.file)), x[[style]] ))
    # pick the species with the highest relative overlap in chromosome names
    gs <- gs[[which.min(sapply(gss, function(x) length(setdiff(x[[2]],x[[1]]))/length(x[[2]])))]]

    s <- gs[[style]][!gs$circular] 
    s[s %in% seqlevels(scanFaIndex(reference.file))]
}    
