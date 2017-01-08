#' Remove low quality targets
#' 
#' This function determines which intervals in the coverage files should be
#' included or excluded in the segmentation. It is called via the
#' \code{fun.filterTargets} argument of \code{\link{runAbsoluteCN}}. The
#' arguments are passed via \code{args.filterTargets}.
#' 
#' 
#' @param log.ratio Copy number log-ratios, one for each target or interval in
#' coverage file.
#' @param tumor GATK coverage data for tumor sample.
#' @param gc.data \code{data.frame} with GC bias for each interval.
#' @param seg.file If not \code{NULL}, then do not filter targets, because data
#' is already segmented via the provided segmentation file.
#' @param filter.lowhigh.gc Quantile q (defines lower q and upper 1-q) for
#' removing targets with outlier GC profile. Assuming that GC correction might
#' not have been worked on those. Requires \code{gc.gene.file}.
#' @param min.targeted.base Exclude intervals with targeted base (size in bp)
#' smaller than this cutoff. This is useful when the same interval file was
#' used to calculate GC content. For such small targets, the GC content is
#' likely very different from the true GC content of the probes.
#' @param normalDB Normal database, created with
#' \code{\link{createNormalDatabase}}.
#' @param normalDB.min.coverage Exclude targets with coverage lower than 20
#' percent of the chromosome median in the pool of normals.
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
#' gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
#'     package="PureCN")
#' 
#' # The max.candidate.solutions, max.ploidy and test.purity parameters are set to
#' # non-default values to speed-up this example.  This is not a good idea for real
#' # samples.
#' ret <-runAbsoluteCN(normal.coverage.file=normal.coverage.file,
#'     tumor.coverage.file=tumor.coverage.file, genome="hg19", vcf.file=vcf.file,
#'     sampleid="Sample1", gc.gene.file=gc.gene.file, normalDB=normalDB,
#'     args.filterTargets=list(min.targeted.base=10), max.ploidy=4, 
#'     test.purity=seq(0.3,0.7,by=0.05), max.candidate.solutions=1)
#' 
#' @export filterTargets
filterTargets <- function(log.ratio, tumor, gc.data, seg.file, 
    filter.lowhigh.gc = 0.001, min.targeted.base = 5, normalDB = NULL,
    normalDB.min.coverage = 0.2) {
    # NA's in log.ratio confuse the CBS function
    targetsUsed <- !is.na(log.ratio) & !is.infinite(log.ratio) 
    # With segmentation file, ignore all filters
    if (!is.null(seg.file)) return(targetsUsed)

    if (!is.null(gc.data)) {
        .checkFraction(filter.lowhigh.gc, "filter.lowhigh.gc")
        targetsUsed <- .filterTargetsLowHighGC(targetsUsed, tumor,
            gc.data, filter.lowhigh.gc)
    }
    targetsUsed <- .filterTargetsNormalDB(targetsUsed, tumor, normalDB,
        normalDB.min.coverage)
    targetsUsed <- .filterTargetsTargetedBase(targetsUsed, tumor,
        min.targeted.base)
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
    tmp <- readCoverageGatk(normalDB$normal.coverage.files[1])
    return(identical(tmp$Target, tumor$Target))
}
    
.filterTargetsNormalDB <- function(targetsUsed, tumor, normalDB,
normalDB.min.coverage) {
    if (is.null(normalDB)) return(targetsUsed)
    # make sure that normalDB matches tumor
    if (!.checkNormalDB(tumor, normalDB)) {
        warning("normalDB does not align with coverage. Ignoring normalDB.")
        return(targetsUsed)
    }    
# TODO: remove in 1.6
    if (!is.null(normalDB$gatk.normal.files)) {
        normalDB$normal.coverage.files <- normalDB$gatk.normal.files
    }    
    nBefore <- sum(targetsUsed)
    min.coverage <- (sapply(split(normalDB$exon.median.coverage, 
        tumor$chr), median, na.rm=TRUE)*normalDB.min.coverage)[tumor$chr]
    
    # should not happen, but just in case
    min.coverage[is.na(min.coverage)] <- median(min.coverage, na.rm=TRUE)
    
    targetsUsed <- targetsUsed & !is.na(normalDB$exon.median.coverage) & 
        normalDB$exon.median.coverage >= min.coverage

    nAfter <- sum(targetsUsed)

    if (nAfter < nBefore) { 
        flog.info("Removing %i targets with low coverage in normalDB.", 
            nBefore-nAfter)
    }
    targetsUsed
}

.filterTargetsChrHash <- function(targetsUsed, tumor, chr.hash) {
    if (is.null(chr.hash)) return(targetsUsed)
    nBefore <- sum(targetsUsed)
    targetsUsed <- targetsUsed & tumor$chr %in% chr.hash$chr
    nAfter <- sum(targetsUsed)

    if ( nAfter < nBefore) { 
        flog.info("Removing %i targets on chromosomes outside chr.hash.", 
            nBefore-nAfter)
    }
    targetsUsed
}

.filterTargetsTargetedBase <- function(targetsUsed, tumor, min.targeted.base) {
    if (is.null(min.targeted.base)) return(targetsUsed)
    nBefore <- sum(targetsUsed)
    targetsUsed <- targetsUsed & !is.na(tumor$targeted.base) & 
        tumor$targeted.base >= min.targeted.base
    nAfter <- sum(targetsUsed)
    if (nAfter < nBefore) {
        flog.info("Removing %i small targets.", nBefore-nAfter)
    }    
    targetsUsed
}

.filterTargetsLowHighGC <- function(targetsUsed, tumor, gc.data,
    filter.lowhigh.gc) {
    gc.data <- gc.data[match(as.character(tumor[,1]), gc.data[,1]),]
    qq <- quantile(gc.data$gc_bias, p=c(filter.lowhigh.gc, 
        1-filter.lowhigh.gc), na.rm=TRUE)

    nBefore <- sum(targetsUsed)
    targetsUsed <- targetsUsed & !is.na(gc.data$gc_bias) & 
        !(gc.data$gc_bias < qq[1] | gc.data$gc_bias > qq[2])
    nAfter <- sum(targetsUsed)

    if (nAfter < nBefore) {
        flog.info("Removing %i low/high GC targets.", nBefore-nAfter)
    }
    targetsUsed
}
