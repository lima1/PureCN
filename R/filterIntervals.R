#' Remove low quality intervals
#' 
#' This function determines which intervals in the coverage files should be
#' included or excluded in the segmentation. It is called via the
#' \code{fun.filterIntervals} argument of \code{\link{runAbsoluteCN}}. The
#' arguments are passed via \code{args.filterIntervals}.
#' 
#' 
#' @param normal Coverage data for normal sample.
#' @param tumor Coverage data for tumor sample.
#' @param log.ratio Copy number log-ratios, one for each interval in the
#' coverage file.
#' @param seg.file If not \code{NULL}, then do not filter intervals, because data
#' is already segmented via the provided segmentation file.
#' @param filter.lowhigh.gc Quantile q (defines lower q and upper 1-q) for
#' removing intervals with outlier GC profile. Assuming that GC correction might
#' not have been worked on those. Requires \code{interval.file}.
#' @param min.coverage Minimum coverage in both normal and tumor. Intervals with
#' lower coverage are ignored. If a \code{normalDB} is provided, then this
#' database already provides information about low quality intervals and the
#' \code{min.coverage} is set to \code{min.coverage/10000}.
#' @param min.total.counts Exclude intervals with fewer than that many reads
#' in combined tumor and normal. 
#' @param min.targeted.base Exclude intervals with targeted base (size in bp)
#' smaller than this cutoff. This is useful when the same interval file was
#' used to calculate GC content. For such small targets, the GC content is
#' likely very different from the true GC content of the probes.
#' @param min.mappability \code{double(2)} specifying the minimum mappability score
#' for on-target, off-target in that order.
#' @param min.fraction.offtarget Skip off-target regions when less than the 
#' specified fraction of all intervals passes all filters
#' @param normalDB Normal database, created with
#' \code{\link{createNormalDatabase}}.
#' @return \code{logical(length(log.ratio))} specifying which intervals should be
#' used in segmentation.
#' @author Markus Riester
#' @examples
#' 
#' normal.coverage.file <- system.file("extdata", "example_normal.txt", 
#'     package="PureCN")
#' normal2.coverage.file <- system.file("extdata", "example_normal2.txt", 
#'     package="PureCN")
#' normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)
#' normalDB <- createNormalDatabase(normal.coverage.files)
#' 
#' tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
#'     package="PureCN")
#' vcf.file <- system.file("extdata", "example.vcf.gz", 
#'     package="PureCN")
#' interval.file <- system.file("extdata", "example_intervals.txt", 
#'     package="PureCN")
#' 
#' # The max.candidate.solutions, max.ploidy and test.purity parameters are set to
#' # non-default values to speed-up this example.  This is not a good idea for real
#' # samples.
#' ret <-runAbsoluteCN(normal.coverage.file=normal.coverage.file,
#'     tumor.coverage.file=tumor.coverage.file, genome="hg19", vcf.file=vcf.file,
#'     sampleid="Sample1", interval.file=interval.file, normalDB=normalDB,
#'     args.filterIntervals=list(min.targeted.base=10), max.ploidy=4, 
#'     test.purity=seq(0.3,0.7,by=0.05), max.candidate.solutions=1)
#' 
#' @export filterIntervals
filterIntervals <- function(normal, tumor, log.ratio, seg.file, 
    filter.lowhigh.gc = 0.001, min.coverage = 15, min.total.counts = 100,
    min.targeted.base = 5, min.mappability = c(0.6, 0.1), 
    min.fraction.offtarget = 0.05, normalDB = NULL) {
    intervalsUsed <- .filterIntervalsNotNA(log.ratio)
    # With segmentation file, ignore all filters
    if (!is.null(seg.file)) return(intervalsUsed)

    if (!is.null(tumor$gc_bias)) {
        .checkFraction(filter.lowhigh.gc, "filter.lowhigh.gc")
        intervalsUsed <- .filterIntervalsLowHighGC(intervalsUsed, tumor,
            filter.lowhigh.gc)
    }
    intervalsUsed <- .filterIntervalsNormalDB(intervalsUsed, tumor, normalDB)

    intervalsUsed <- .filterIntervalsTargetedBase(intervalsUsed, tumor,
        min.targeted.base)
    intervalsUsed <- .filterIntervalsTotalNormalCoverage(intervalsUsed, normal,
        min.targeted.base, min.coverage)

    if (!is.null(normalDB)) {
        min.coverage <- min.coverage/10000
        flog.info("normalDB provided. Setting minimum coverage for segmentation to %.4fX.", min.coverage)
    } else {
        flog.warn("No normalDB provided. Provide one for better results.")
    }    
        
    intervalsUsed <- .filterIntervalsCoverage(intervalsUsed, normal, tumor, 
        min.coverage)

    intervalsUsed <- .filterIntervalsTotalCounts(intervalsUsed, normal, tumor, 
        min.total.counts)

    intervalsUsed <- .filterIntervalsMappability(intervalsUsed, tumor,
        min.mappability)
    
    intervalsUsed <- .filterIntervalsOfftarget(intervalsUsed, tumor,
                                               min.fraction.offtarget)
    return(intervalsUsed)
}

    
.checkNormalDB <- function(tumor, normalDB) {
    if (!is(normalDB, "list")) {
        .stopUserError("normalDB not a valid normalDB object. ",
            "Use createNormalDatabase to create one.")
    }    
    if (is.null(normalDB$version) || normalDB$version < 4) {
        .stopUserError("normalDB incompatible with this PureCN version. ",
                       "Please re-run NormalDB.R.")
    }
    if (is.null(normalDB$normal.coverage.files) ||
        !length(normalDB$normal.coverage.files)) {
        .stopUserError("normalDB appears to be empty.")
    }    
    intervals <- normalDB[["intervals"]]
    # TODO remove in PureCN 1.14 
    if (is.null(intervals)) {
        intervals <- as.character(readCoverageFile(normalDB$normal.coverage.files[1]))
    }
    return(identical(intervals, as.character(tumor)))
}

.filterIntervalsNotNA <- function(log.ratio) {
    nBefore <- length(log.ratio) 
    # NA's in log.ratio confuse the CBS function
    intervalsUsed <- !is.na(log.ratio) & !is.infinite(log.ratio) 
    nAfter <- sum(intervalsUsed)

    if (nAfter < nBefore) {
        flog.info("Removing %i intervals with missing log.ratio.", 
            nBefore-nAfter)
    }

    intervalsUsed
}
        
.filterIntervalsNormalDB <- function(intervalsUsed, tumor, normalDB,
normalDB.min.coverage, normalDB.max.missing) {
    if (is.null(normalDB)) return(intervalsUsed)
    # make sure that normalDB matches tumor
    if (!.checkNormalDB(tumor, normalDB)) {
        flog.warn("normalDB does not align with coverage. Ignoring normalDB.")
        return(intervalsUsed)
    }    

    nBefore <- sum(intervalsUsed) 
    intervalsUsed <- intervalsUsed & normalDB$intervals.used
    nAfter <- sum(intervalsUsed)

    if (nAfter < nBefore) {
        flog.info("Removing %i intervals excluded in normalDB.", 
            nBefore-nAfter)
    }
    intervalsUsed
}

.filterIntervalsChrHash <- function(intervalsUsed, tumor, chr.hash) {
    if (is.null(chr.hash)) return(intervalsUsed)
    nBefore <- sum(intervalsUsed)
    intervalsUsed <- intervalsUsed & seqnames(tumor) %in% chr.hash$chr
    nAfter <- sum(intervalsUsed)

    if (nAfter < nBefore) {
        flog.info("Removing %i intervals on chromosomes outside chr.hash.", 
            nBefore-nAfter)
    }
    intervalsUsed
}

.filterIntervalsTargetedBase <- function(intervalsUsed, tumor, min.targeted.base) {
    if (is.null(min.targeted.base)) return(intervalsUsed)
    nBefore <- sum(intervalsUsed)
    intervalsUsed <- intervalsUsed & width(tumor) >= min.targeted.base
    nAfter <- sum(intervalsUsed)
    if (nAfter < nBefore) {
        flog.info("Removing %i small (< %ibp) intervals.", nBefore-nAfter, 
            min.targeted.base)
    }    
    intervalsUsed
}

.filterIntervalsTotalNormalCoverage <- function(intervalsUsed, normal,
    min.targeted.base, min.coverage) {
    if (is.null(min.targeted.base) || is.null(min.coverage)) {
        return(intervalsUsed)
    }    
    nBefore <- sum(intervalsUsed)
    cutoff <- min.targeted.base * min.coverage * 2
    intervalsUsed <- intervalsUsed & normal$coverage >= cutoff
    nAfter <- sum(intervalsUsed)
    if (nAfter < nBefore) {
        flog.info("Removing %i intervals with low total coverage in normal (< %.2f reads).", 
            nBefore-nAfter, cutoff)
    }    
    intervalsUsed
}

.filterIntervalsCoverage <- function(intervalsUsed, normal, tumor, min.coverage) {
    #MR: we try to not remove homozygous deletions in very pure samples.
    #  to distinguish low quality from low copy number, we keep if normal
    # has good coverage. If normal coverage is very high, we adjust for that.  
    total.cov.normal <- sum(as.numeric(normal$coverage), na.rm = TRUE)
    total.cov.tumor <- sum(as.numeric(tumor$coverage), na.rm = TRUE)

    f <- max(total.cov.normal / total.cov.tumor, 1)

    well.covered.interval.idx <- (normal$average.coverage >= min.coverage &
        tumor$average.coverage >= min.coverage) | 
        (normal$average.coverage >= 1.5 * f * min.coverage &  
        tumor$average.coverage >= 0.5 * min.coverage)
    #MR: fix for missing chrX/Y
    well.covered.interval.idx[is.na(well.covered.interval.idx)] <- FALSE

    nBefore <- sum(intervalsUsed)
    intervalsUsed <- intervalsUsed & well.covered.interval.idx
    nAfter <- sum(intervalsUsed)
    if (nAfter < nBefore) {
        flog.info("Removing %i low coverage (< %.4fX) intervals.", nBefore-nAfter, 
            min.coverage)
    }    
    intervalsUsed
}
.filterIntervalsTotalCounts <- function(intervalsUsed, normal, tumor, min.total.counts) {

    well.covered.interval.idx <- ( normal$counts + tumor$counts ) >= min.total.counts
    well.covered.interval.idx[is.na(well.covered.interval.idx)] <- FALSE

    nBefore <- sum(intervalsUsed)
    intervalsUsed <- intervalsUsed & well.covered.interval.idx
    nAfter <- sum(intervalsUsed)
    if (nAfter < nBefore) {
        flog.info("Removing %i low count (< %.0f total reads) intervals.", nBefore-nAfter, 
            min.total.counts)
    }    
    intervalsUsed
}


.filterIntervalsLowHighGC <- function(intervalsUsed, tumor, filter.lowhigh.gc) {
    qq <- quantile(tumor$gc_bias, p=c(filter.lowhigh.gc, 
        1-filter.lowhigh.gc), na.rm=TRUE)

    nBefore <- sum(intervalsUsed)
    intervalsUsed <- intervalsUsed & !is.na(tumor$gc_bias) & 
        !(tumor$gc_bias < qq[1] | tumor$gc_bias > qq[2])
    nAfter <- sum(intervalsUsed)

    if (nAfter < nBefore) {
        flog.info("Removing %i low/high GC targets.", nBefore-nAfter)
    }
    intervalsUsed
}

.filterIntervalsMappability <- function(intervalsUsed, tumor, min.mappability) {
    if (is.null(tumor$mappability) || all(is.na(tumor$mappability))) {
        return(intervalsUsed)
    }    
    if (any(is.na(tumor$mappability))) {
        flog.warn("Some intervals do not have a mappability score, assuming they are fine.")
    }
    nBefore <- sum(intervalsUsed)
    intervalsUsed <- intervalsUsed & ( is.na(tumor$mappability) |
        (tumor$on.target & tumor$mappability >= min.mappability[1]) |
        (!tumor$on.target & tumor$mappability >= min.mappability[2]))
    nAfter <- sum(intervalsUsed)

    if (nAfter < nBefore) {
        flog.info("Removing %i low mappability intervals.", nBefore-nAfter)
    }
    intervalsUsed
}

.filterIntervalsOfftarget <- function(intervalsUsed, tumor,
                                      min.fraction.offtarget) {
    
    n <- sum(!tumor$on.target[intervalsUsed], na.rm = TRUE)
    m <- sum(tumor$on.target[intervalsUsed], na.rm = TRUE)
    f <- n / (n + m)
    if (!n) {
        return(intervalsUsed)
    }    
    nBefore <- sum(intervalsUsed)
    if (f < min.fraction.offtarget) {
        flog.warn("Not enough off-target intervals. Ignoring them (%i on-target, %i off-target, %.2f).",
            n, m, f)
        intervalsUsed <- intervalsUsed & !is.na(tumor$on.target) &
            tumor$on.target
    } else if (f < 0.1) {
        flog.warn("Low number of off-target intervals. You might want to exclude them (%i on-target, %i off-target, %.2f).",
            n, m, f)
    }    
    nAfter <- sum(intervalsUsed)

    if (nAfter < nBefore) {
        flog.info("Removing %i off-target intervals.", nBefore-nAfter)
    }
    intervalsUsed
}
