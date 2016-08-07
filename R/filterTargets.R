filterTargets <- structure(function(# Remove low quality targets
### This function determines which intervals in the coverage files should
### be included or excluded in the segmentation. It is called via the
### \code{fun.filterTargets} argument of \code{\link{runAbsoluteCN}}. The 
### arguments are passed via \code{args.filterTargets}.
log.ratio, 
### Copy number log-ratios, one for each target or interval
### in coverage file.
tumor, 
### GATK coverage data for tumor sample.
gc.data, 
### \code{data.frame} with GC bias for each interval.
seg.file,
### If not \code{NULL}, then do not filter targets, because data is
### already segmented via the provided segmentation file.
filter.lowhigh.gc=0.001,
### Quantile q (defines lower q and upper 1-q) 
### for removing targets with outlier GC profile. Assuming that GC correction 
### might not have been worked on those. Requires \code{gc.gene.file}.
filter.targeted.base=4,
### Exclude exons with targeted base (size) smaller 
### than this cutoff. This is useful when the same interval file was used to
### calculate GC content. For such small exons, the GC content is likely 
### very different from the true GC content of the probes.
normalDB=NULL,
### Normal database, created with \code{\link{createNormalDatabase}}.
normalDB.min.coverage=0.2,
### Exclude targets with coverage lower than 20 percent of the 
### chromosome median in the pool of normals.
verbose
### Verbose output.
) {
    # NA's in log.ratio confuse the CBS function
    targetsUsed <- !is.na(log.ratio) & !is.infinite(log.ratio) 
    # With segmentation file, ignore all filters
    if (!is.null(seg.file)) return(targetsUsed)

    if (!is.null(gc.data)) {
        .checkFraction(filter.lowhigh.gc, "filter.lowhigh.gc")
        targetsUsed <- .filterTargetsLowHighGC(targetsUsed, tumor,
            gc.data, filter.lowhigh.gc, verbose)
    }
    targetsUsed <- .filterTargetsNormalDB(targetsUsed, tumor, normalDB,
        normalDB.min.coverage, verbose)
    targetsUsed <- .filterTargetsTargetedBase(targetsUsed, tumor,
        filter.targeted.base, verbose)
},ex=function() {
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gatk.normal2.file <- system.file("extdata", "example_normal2.txt", 
    package="PureCN")
gatk.normal.files <- c(gatk.normal.file, gatk.normal2.file)
normalDB <- createNormalDatabase(gatk.normal.files)

gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")
vcf.file <- system.file("extdata", "example_vcf.vcf", 
    package="PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
    package="PureCN")

# The max.candidate.solutions, max.ploidy and test.purity parameters are set to
# non-default values to speed-up this example.  This is not a good idea for real
# samples.
ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file,
    gatk.tumor.file=gatk.tumor.file, genome="hg19", vcf.file=vcf.file,
    sampleid='Sample1', gc.gene.file=gc.gene.file,
    args.filterTargets=list(normalDB=normalDB), max.ploidy=4, 
    test.purity=seq(0.3,0.7,by=0.05), max.candidate.solutions=1)
})

.filterTargetsNormalDB <- function(targetsUsed, tumor, normalDB,
normalDB.min.coverage, verbose) {
    if (is.null(normalDB)) return(targetsUsed)
    nBefore <- sum(targetsUsed)
    min.coverage <- (sapply(split(normalDB$exon.median.coverage, 
        tumor$chr), median, na.rm=TRUE)*normalDB.min.coverage)[tumor$chr]
    targetsUsed <- targetsUsed & tumor$average.coverage >= min.coverage
    nAfter <- sum(targetsUsed)

    if (verbose && nAfter < nBefore) { 
        message("Removing ", nBefore-nAfter, " exons with low coverage ",
            "in normalDB.")
    }
    targetsUsed
}

.filterTargetsChrHash <- function(targetsUsed, tumor, chr.hash, verbose) {
    if (is.null(chr.hash)) return(targetsUsed)
    nBefore <- sum(targetsUsed)
    targetsUsed <- targetsUsed & tumor$chr %in% chr.hash$chr
    nAfter <- sum(targetsUsed)

    if (verbose && nAfter < nBefore) { 
        message("Removing ", nBefore-nAfter, " exons on chromosomes ",
            "outside chr.hash.")
    }
    targetsUsed
}
.filterTargetsTargetedBase <- function(targetsUsed, tumor, filter.targeted.base,
    verbose) {
    if (is.null(filter.targeted.base)) return(targetsUsed)
    nBefore <- sum(targetsUsed)
    targetsUsed <- targetsUsed & !is.na(tumor$targeted.base) & 
        tumor$targeted.base >= filter.targeted.base
    nAfter <- sum(targetsUsed)
    if (verbose && nAfter < nBefore) message("Removing ", nBefore-nAfter, 
        " small exons.")
    targetsUsed
}
.filterTargetsLowHighGC <- function(targetsUsed, tumor, gc.data,
    filter.lowhigh.gc, verbose) {
    gc.data <- gc.data[match(as.character(tumor[,1]), gc.data[,1]),]
    qq <- quantile(gc.data$gc_bias, p=c(filter.lowhigh.gc, 
        1-filter.lowhigh.gc), na.rm=TRUE)

    nBefore <- sum(targetsUsed)
    targetsUsed <- targetsUsed & !is.na(gc.data$gc_bias) & 
        !(gc.data$gc_bias < qq[1] | gc.data$gc_bias > qq[2])
    nAfter <- sum(targetsUsed)

    if (verbose && nAfter < nBefore) message("Removing ", 
        nBefore-nAfter, " low/high GC exons.")

    targetsUsed
}
