#' Preprocess intervals
#' 
#' Optimize intervals for copy number calling by tiling long intervals and by 
#' including off-target regions. Uses \code{scanFa} from the Rsamtools package 
#' to retrieve GC content of intervals in a reference FASTA file. If provided,
#' will annotate intervals with mappability and replication timing scores.
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
#' @param min.target.width Make sure that target regions are of at least
#' this specified width. See \code{small.targets}.
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
#' @param reptiming Annotate intervals with replication timing score. Expected as 
#' \code{GRanges} object with first meta column being the score. 
#' @param average.reptiming.width Tile \code{reptiming} into bins of specified
#' width. 
#' @param exclude Any target that overlaps with this \code{GRanges} object
#' will be excluded. 
#' @param off.target.seqlevels Controls how to deal with chromosomes/contigs
#' found in the \code{reference.file} but not in the \code{interval.file}.
#' @param small.targets Strategy to deal with targets smaller than
#' \code{min.target.width}.
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
#' preprocessIntervals(interval.file, reference.file, 
#'     output.file="gc_file.txt")
#' 
#' intervals <- import(bed.file)
#' preprocessIntervals(intervals, reference.file, 
#'     output.file="gc_file.txt")
#' 
#' @export preprocessIntervals
#' @importFrom BiocGenerics unstrand score
#' @importFrom Biostrings letterFrequency
#' @importFrom GenomeInfoDb seqlengths seqlevelsInUse seqlevels<- seqlengths<-
#' @importFrom GenomicRanges tileGenome
#' @importFrom S4Vectors mcols
#' @importFrom rtracklayer import
#' @importFrom methods is
#' @importFrom stats aggregate
preprocessIntervals <- function(interval.file, reference.file,
                                output.file = NULL, off.target = FALSE,
                                average.target.width = 400,
                                min.target.width = 100,
                                min.off.target.width = 20000,
                                average.off.target.width = 200000,
                                off.target.padding = -500, mappability = NULL,
                                min.mappability = c(0.5, 0.1, 0.7), 
                                reptiming = NULL,
                                average.reptiming.width = 100000,
                                exclude = NULL,
                                off.target.seqlevels=c("targeted", "all"),
                                small.targets=c("resize", "drop")) {

    if (is(interval.file, "GRanges")) {
        interval.gr <- .checkIntervals(unstrand(interval.file))
    } else {
        interval.gr <- readCoverageFile(interval.file)
    }

    # make sure the chromsome naming style is the same in all provided files
    # be nice and fix it if necessary
    interval.gr <- .checkSeqlevelStyle(scanFaIndex(reference.file), interval.gr, "interval")
    interval.gr <- .checkSeqlengths(scanFaIndex(reference.file), interval.gr)
    interval.gr <- .checkTargetWidth(interval.gr, min.target.width,
        match.arg(small.targets))
    if (is.null(interval.gr$on.target)) interval.gr$on.target <- TRUE    
    if (!is.null(mappability)) {
        mappability <- .checkSeqlevelStyle(scanFaIndex(reference.file), mappability, "mappability")
        mappability <- .remove0MappabilityRegions(mappability)
    }
    if (!is.null(reptiming)) {
        reptiming <- .checkSeqlevelStyle(scanFaIndex(reference.file), reptiming, 
                                         "reptiming")
        reptiming <- .tileReptiming(seqinfo(scanFaIndex(reference.file)), 
                                    reptiming, average.reptiming.width)
    }
     
    containsOfftarget <- sum(interval.gr$on.target)!=length(interval.gr)

    if (containsOfftarget) {
        flog.warn("Intervals contain off-target regions. %s",
            "Will not change intervals.")
    } else {
        interval.gr <- .splitIntervals(interval.gr, average.target.width)
        interval.gr <- .excludeIntervals(interval.gr, exclude)

        # find off-target regions
        if (off.target) {
            off.target.seqlevels <- match.arg(off.target.seqlevels)
            interval.gr <- .annotateIntervalsOfftarget(interval.gr, 
                reference.file, mappability, off.target.padding,
                min.off.target.width, average.off.target.width, 
                off.target.seqlevels)
        }
    }

    interval.gr <- .annotateIntervalsMappability(interval.gr, mappability, 
        min.mappability)
    interval.gr <- .annotateIntervalsReptiming(interval.gr, reptiming)
    interval.gr <- .annotateIntervalsGC(interval.gr, reference.file)
        
    if (is.null(interval.gr$Gene)) interval.gr$Gene <- "."

    if (!is.null(output.file)) {
        .writeIntervals(interval.gr, output.file)
    }    
    invisible(interval.gr)
}

# this function removes short chromosomes that have no probes (mainly a
# general way to remove chrM)
.dropShortUntargetedSeqLevels <- function(offRegions, interval.gr, minSize) {
    idx <- seqlevels(offRegions) %in% seqlevels(interval.gr) |
        seqlengths(offRegions) >= minSize
    offRegions <- offRegions[seqnames(offRegions) %in% seqlevels(offRegions)[idx]]
    seqlevels(offRegions) <- seqlevelsInUse(offRegions)
    offRegions
}

.annotateIntervalsOfftarget <- function(interval.gr, reference.file, 
                                        mappability, off.target.padding,
                                        min.off.target.width,
                                        average.off.target.width, 
                                        off.target.seqlevels) {

    if (off.target.padding > 0) {
        .stopUserError("off.target.padding must be negative.")
    }    
    offRegions <- setdiff(scanFaIndex(reference.file), unstrand(interval.gr))
    offRegions <- .dropShortUntargetedSeqLevels(offRegions, interval.gr, 
        average.off.target.width)

    if (!is.null(mappability)) {
        offRegions <- intersect(offRegions, mappability)
    }    
    offRegions <- offRegions[width(offRegions)>off.target.padding*-2]
    if (!length(offRegions)) {
        .stopUserError("No off-target regions after filtering for mappability ",
            "and off.target.padding")
    }
    offRegions <- .padGranges(offRegions, off.target.padding)

    flog.info("Tiling off-target regions to an average width of %i.",
        average.off.target.width)
        
    offRegions <- unlist(tile(offRegions, width=average.off.target.width))
    offRegions$on.target <- FALSE
    offRegions <- offRegions[width(offRegions)>=min.off.target.width]
    seqlevelsBefore <- seqlevelsInUse(offRegions)
    if (off.target.seqlevels == "targeted") {
        offRegions <- offRegions[seqnames(offRegions) %in% 
            seqlevelsInUse(interval.gr)]
    }
    seqlevelsAfter <- seqlevelsInUse(offRegions)
    if (!identical(seqlevelsBefore, seqlevelsAfter)) {
        flog.info("Removing following contigs from off-target regions: %s", 
            paste(setdiff(seqlevelsBefore, seqlevelsAfter), collapse=","))
    }
    merge(offRegions, interval.gr, all = TRUE, sort=TRUE)
}

# add GC content 
.annotateIntervalsGC <- function(interval.gr, reference.file) {
    flog.info("Calculating GC-content...")
    x <- scanFa(reference.file, interval.gr)
    GC.count <- letterFrequency(x,"GC")
    all.count <- letterFrequency(x,"ATGC")
    interval.gr$gc_bias <- as.vector(ifelse(all.count==0,NA,GC.count/all.count))
    # exclude unavailable regions
    interval.gr[which(!is.na(interval.gr$gc_bias))]
}
    
# add mappability score to intervals
.annotateIntervalsMappability <- function(interval.gr, mappability, min.mappability) {
    interval.gr <- .addScoreToGr(interval.gr, mappability, "mappability")
    if (is.null(mappability)) return(interval.gr)

    # remove intervals with low mappability
    nBefore <- sum(interval.gr$on.target)
    interval.gr <- interval.gr[
        (interval.gr$on.target & interval.gr$mappability >= min.mappability[1]) |
        (!interval.gr$on.target & interval.gr$mappability >= min.mappability[2]) ]
    # remove chrY low mappability    
    sex.chr <- .getSexChr(seqlevels(interval.gr))[2]
    interval.gr <- interval.gr[!seqnames(interval.gr) %in% sex.chr |
        interval.gr$mappability >= min.mappability[3] ]
    nAfter <- sum(interval.gr$on.target)
    if (nBefore > nAfter) {
        flog.info("Removing %i intervals with low mappability score (<%.2f).", 
            nBefore-nAfter, min.mappability[1])
    }
    interval.gr
}

.annotateIntervalsReptiming <- function(interval.gr, reptiming) {
    .addScoreToGr(interval.gr, reptiming, "reptiming")
}

.checkColScore <- function(y, label) {
    colScore <- if (is.null(y$score)) 1 else "score"
    if (!is(mcols(y)[, colScore], "numeric")) {
        flog.warn("Score column in %s file is not numeric.", label)
        class(mcols(y)[, colScore]) <- "numeric"
    }
    y
}
.getColScore <- function(y) {
    colScore <- if (is.null(y$score)) 1 else "score"
}
.addScoreToGr <- function(interval.gr, y, label) {
    mcols(interval.gr)[[label]] <- NA
    if (!is.null(y)) {
        y <- .checkColScore(y, label)
        ov <- findOverlaps(interval.gr, y)
        colScore <- .getColScore(y)

        mappScore <- aggregate(mcols(y)[subjectHits(ov),colScore], 
            by=list(queryHits(ov)), mean)
        mcols(interval.gr)[[label]][mappScore[,1]] <- mappScore[,2]
        idxNA <- is.na(mcols(interval.gr)[[label]])

        if (sum(idxNA)) {
            if (!is.null(interval.gr$on.target)) {
                sumOntarget <- sum(idxNA & interval.gr$on.target, na.rm = TRUE)
                flog.warn("%i intervals without %s score (%i on-target).", 
                    sum(idxNA), label, sumOntarget)
            }
            mcols(interval.gr)[[label]][idxNA] <- 0
        }    
    } else {
        flog.warn("No %s scores provided.", label)
    }    
    return(interval.gr)
}

.splitIntervals <- function(interval.gr, average.target.width) {
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
    interval.gr
}

.excludeIntervals <- function(interval.gr, exclude) {
    if (is.null(exclude)) return(interval.gr)
    nBefore <- length(interval.gr)
    interval.gr <- interval.gr[which(!overlapsAny(interval.gr, exclude) | 
                                     !interval.gr$on.target)]

    nAfter <- length(interval.gr)
    if (nBefore > nAfter) {
        flog.info("Removing %i targets overlapping with exclude.", 
            nBefore - nAfter)
    }
    interval.gr
}

.remove0MappabilityRegions <- function(mappability) {
    mappability <- .checkColScore(mappability, "mappability")
    colScore <- .getColScore(mappability)
    mappability[which(mcols(mappability)[, colScore]>0),]
}
        
.writeIntervals <- function(interval.gr, output.file) {
    tmp <- data.frame(
        Target=as.character(interval.gr),
        gc_bias=interval.gr$gc_bias,
        mappability=interval.gr$mappability,
        reptiming=interval.gr$reptiming,
        Gene=interval.gr$Gene,
        on_target=interval.gr$on.target
    )    
    fwrite(tmp, file = output.file, row.names = FALSE, quote = FALSE, sep = "\t")
}

.checkSeqlengths <- function(ref, x) {
    isc <- intersect(names(seqlengths(x)), names(seqlengths(ref)))
    if (length(isc)) {
        seqlengths(x)[isc] <- seqlengths(ref)[isc]
    }
    x
}    
.checkSeqlevelStyle <- function(ref, x, name1, name2="reference") {
    refSeqlevelStyle <- try(.getSeqlevelsStyle(ref), silent=TRUE)
    # if unknown, we cannot check and correct
    if (is(refSeqlevelStyle, "try-error")) return(x)
    xSeqlevelStyle <- try(.getSeqlevelsStyle(x), silent=TRUE)

    if (is(xSeqlevelStyle, "try-error")) {
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
.tileReptiming <- function(seqinfo.ref, reptiming, average.reptiming.width) {
    if (is.null(average.reptiming.width) || average.reptiming.width <= 1) {
        return(reptiming)
    }
    flog.info("Averaging reptiming into bins of size %i...", average.reptiming.width)
    bins <- tileGenome(seqinfo.ref, tilewidth = average.reptiming.width, 
        cut.last.tile.in.chrom = TRUE)
    reptiming <- .addScoreToGr(bins, reptiming, "score")
    reptiming <- reptiming[!is.na(score(reptiming))]
}

.checkTargetWidth <- function(interval.gr, min.target.width, small.targets) {
    idx <- width(interval.gr) < min.target.width
    if (any(idx)) {
        flog.warn("Found small target regions (< %ibp). Will %s them.",
            min.target.width, small.targets)
        if (small.targets == "drop") {
            interval.gr <- interval.gr[!idx,]
        } else {
            off <- floor((width(interval.gr[idx]) - min.target.width) / 2)
            start(interval.gr[idx]) <- pmax(start(interval.gr[idx]) + off, 1)
            end(interval.gr[idx]) <- pmin(start(interval.gr[idx]) + min.target.width - 1,
                seqlengths(interval.gr)[as.character(seqnames(interval.gr[idx]))], 
                na.rm = TRUE)
            interval.gr <- reduce(interval.gr)
        }
    }
    interval.gr
}    
