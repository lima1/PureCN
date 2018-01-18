#' Remove low quality targets
#' 
#' This function determines which intervals in the coverage files should be
#' included or excluded in the segmentation. It is called via the
#' \code{fun.filterTargets} argument of \code{\link{runAbsoluteCN}}. The
#' arguments are passed via \code{args.filterTargets}.
#' 
#' 
#' @param normal Coverage data for normal sample.
#' @param tumor Coverage data for tumor sample.
#' @param log.ratio Copy number log-ratios, one for each target or interval in
#' coverage file.
#' @param seg.file If not \code{NULL}, then do not filter targets, because data
#' is already segmented via the provided segmentation file.
#' @param filter.lowhigh.gc Quantile q (defines lower q and upper 1-q) for
#' removing targets with outlier GC profile. Assuming that GC correction might
#' not have been worked on those. Requires \code{interval.file}.
#' @param min.coverage Minimum coverage in both normal and tumor. Targets with
#' lower coverage are ignored. If a \code{normalDB} is provided, then this
#' database already provides information about low quality targets and the
#' \code{min.coverage} is set to \code{min.coverage/10000}.
#' @param min.targeted.base Exclude intervals with targeted base (size in bp)
#' smaller than this cutoff. This is useful when the same interval file was
#' used to calculate GC content. For such small targets, the GC content is
#' likely very different from the true GC content of the probes.
#' @param normalDB Normal database, created with
#' \code{\link{createNormalDatabase}}.
#' @return \code{logical(length(log.ratio))} specifying which targets should be
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
#' vcf.file <- system.file("extdata", "example_vcf.vcf", 
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
#'     args.filterTargets=list(min.targeted.base=10), max.ploidy=4, 
#'     test.purity=seq(0.3,0.7,by=0.05), max.candidate.solutions=1)
#' 
#' @export filterTargets
filterTargets <- function(normal, tumor, log.ratio, seg.file, 
    filter.lowhigh.gc = 0.001, min.coverage = 15, min.targeted.base = 5, 
    normalDB = NULL) {
    targetsUsed <- .filterTargetsNotNA(log.ratio)
    # With segmentation file, ignore all filters
    if (!is.null(seg.file)) return(targetsUsed)

    if (!is.null(tumor$gc_bias)) {
        .checkFraction(filter.lowhigh.gc, "filter.lowhigh.gc")
        targetsUsed <- .filterTargetsLowHighGC(targetsUsed, tumor,
            filter.lowhigh.gc)
    }
    targetsUsed <- .filterTargetsNormalDB(targetsUsed, tumor, normalDB)

    targetsUsed <- .filterTargetsTargetedBase(targetsUsed, tumor,
        min.targeted.base)
    targetsUsed <- .filterTargetsTotalNormalCoverage(targetsUsed, normal,
        min.targeted.base, min.coverage)

    if (!is.null(normalDB)) {
        min.coverage <- min.coverage/10000
        flog.info("normalDB provided. Setting minimum coverage for segmentation to %.4fX.", min.coverage)
    } else {
        flog.warn("No normalDB provided. Provide one for better results.")
    }    
        
    targetsUsed <- .filterTargetsCoverage(targetsUsed, normal, tumor, 
        min.coverage)

    return(targetsUsed)
}


.checkNormalDB <- function(tumor, normalDB) {
    if (!class(normalDB) == "list") {
        .stopUserError("normalDB not a valid normalDB object. ",
            "Use createNormalDatabase to create one.")
    }    
# TODO: remove in 1.6
    if (!is.null(normalDB$gatk.normal.files)) {
        normalDB$normal.coverage.files <- normalDB$gatk.normal.files
    }    
    if (is.null(normalDB$normal.coverage.files) ||
        !length(normalDB$normal.coverage.files)) {
        .stopUserError("normalDB appears to be empty.")
    }    
    tmp <- readCoverageFile(normalDB$normal.coverage.files[1])
    return(identical(tmp$Target, tumor$Target))
}

.filterTargetsNotNA <- function(log.ratio) {
    nBefore <- length(log.ratio) 
    # NA's in log.ratio confuse the CBS function
    targetsUsed <- !is.na(log.ratio) & !is.infinite(log.ratio) 
    nAfter <- sum(targetsUsed)

    if (nAfter < nBefore) {
        flog.info("Removing %i intervals with missing log.ratio.", 
            nBefore-nAfter)
    }

    targetsUsed
}
        
.filterTargetsNormalDB <- function(targetsUsed, tumor, normalDB,
normalDB.min.coverage, normalDB.max.missing) {
    if (is.null(normalDB)) return(targetsUsed)
    # make sure that normalDB matches tumor
    if (!.checkNormalDB(tumor, normalDB)) {
        flog.warn("normalDB does not align with coverage. Ignoring normalDB.")
        return(targetsUsed)
    }    

    nBefore <- sum(targetsUsed) 
    targetsUsed <- targetsUsed & normalDB$intervals.used
    nAfter <- sum(targetsUsed)

    if (nAfter < nBefore) {
        flog.info("Removing %i targets excluded in normalDB.", 
            nBefore-nAfter)
    }
    targetsUsed
}

.filterTargetsChrHash <- function(targetsUsed, tumor, chr.hash) {
    if (is.null(chr.hash)) return(targetsUsed)
    nBefore <- sum(targetsUsed)
    targetsUsed <- targetsUsed & seqnames(tumor) %in% chr.hash$chr
    nAfter <- sum(targetsUsed)

    if (nAfter < nBefore) {
        flog.info("Removing %i targets on chromosomes outside chr.hash.", 
            nBefore-nAfter)
    }
    targetsUsed
}

.filterTargetsTargetedBase <- function(targetsUsed, tumor, min.targeted.base) {
    if (is.null(min.targeted.base)) return(targetsUsed)
    nBefore <- sum(targetsUsed)
    targetsUsed <- targetsUsed & width(tumor) >= min.targeted.base
    nAfter <- sum(targetsUsed)
    if (nAfter < nBefore) {
        flog.info("Removing %i small (< %ibp) targets.", nBefore-nAfter, 
            min.targeted.base)
    }    
    targetsUsed
}

.filterTargetsTotalNormalCoverage <- function(targetsUsed, normal,
    min.targeted.base, min.coverage) {
    if (is.null(min.targeted.base) || is.null(min.coverage)) {
        return(targetsUsed)
    }    
    nBefore <- sum(targetsUsed)
    cutoff <- min.targeted.base * min.coverage * 2
    targetsUsed <- targetsUsed & normal$coverage >= cutoff
    nAfter <- sum(targetsUsed)
    if (nAfter < nBefore) {
        flog.info("Removing %i targets with low total coverage in normal (< %.2f reads).", 
            nBefore-nAfter, cutoff)
    }    
    targetsUsed
}

.filterTargetsCoverage <- function(targetsUsed, normal, tumor, min.coverage) {
    #MR: we try to not remove homozygous deletions in very pure samples.
    #  to distinguish low quality from low copy number, we keep if normal
    # has good coverage. If normal coverage is very high, we adjust for that.  
    total.cov.normal <- sum(as.numeric(normal$coverage), na.rm = TRUE)
    total.cov.tumor <- sum(as.numeric(tumor$coverage), na.rm = TRUE)

    f <- max(total.cov.normal / total.cov.tumor, 1)

    well.covered.exon.idx <- (normal$average.coverage >= min.coverage &
        tumor$average.coverage >= min.coverage) | 
        (normal$average.coverage >= 1.5 * f * min.coverage &  
        tumor$average.coverage >= 0.5 * min.coverage)
    #MR: fix for missing chrX/Y
    well.covered.exon.idx[is.na(well.covered.exon.idx)] <- FALSE

    nBefore <- sum(targetsUsed)
    targetsUsed <- targetsUsed & well.covered.exon.idx
    nAfter <- sum(targetsUsed)
    if (nAfter < nBefore) {
        flog.info("Removing %i low coverage (< %.4fX) targets.", nBefore-nAfter, 
            min.coverage)
    }    
    targetsUsed
}

.filterTargetsLowHighGC <- function(targetsUsed, tumor, filter.lowhigh.gc) {
    qq <- quantile(tumor$gc_bias, p=c(filter.lowhigh.gc, 
        1-filter.lowhigh.gc), na.rm=TRUE)

    nBefore <- sum(targetsUsed)
    targetsUsed <- targetsUsed & !is.na(tumor$gc_bias) & 
        !(tumor$gc_bias < qq[1] | tumor$gc_bias > qq[2])
    nAfter <- sum(targetsUsed)

    if (nAfter < nBefore) {
        flog.info("Removing %i low/high GC targets.", nBefore-nAfter)
    }
    targetsUsed
}
