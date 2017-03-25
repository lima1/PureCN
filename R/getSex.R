#' Get sample sex from coverage
#' 
#' This function determines the sex of a sample by the coverage ratio of chrX
#' and chrY. Loss of chromosome Y (LOY) can result in a wrong female call. For
#' small targeted panels, this will only work when sufficient sex marker genes
#' such as AMELY are covered. For optimal results, parameters might need to be
#' tuned for the assay.
#' 
#' 
#' @param coverage.file GATK coverage file or data read with
#' \code{\link{readCoverageFile}}.
#' @param min.ratio Min chrX/chrY coverage ratio to call sample as female.
#' @param min.ratio.na Min chrX/chrY coverage ratio to call sample as
#' \code{NA}.  This ratio defines a grey zone from \code{min.ratio.na} to
#' \code{min.ratio} in which samples are not called. The default is set to a
#' copy number ratio that would be rare in male samples, but lower than
#' expected in female samples. Contamination can be a source of ambiguous
#' calls. Mappability issues on chromosome Y resulting in low coverage need to
#' be considered when setting cutoffs.
#' @param remove.outliers Removes coverage outliers before calculating mean
#' chromosome coverages.
#' @return Returns a \code{character(1)} with \code{M} for male, \code{F} for
#' female, or \code{NA} if unknown.
#' @author Markus Riester
#' @seealso \code{\link{getSexFromVcf}}
#' @examples
#' 
#' tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
#'     package="PureCN")
#' sex <- getSexFromCoverage(tumor.coverage.file)
#' 
#' @export getSexFromCoverage
getSexFromCoverage <- function(coverage.file, min.ratio = 25, min.ratio.na = 20,
    remove.outliers = TRUE) {
    if (is.character(coverage.file)) {
        x <- readCoverageFile(coverage.file)
    } else {
        x <- coverage.file
    }

    sex.chr <- .getSexChr(x$chr)
    xx <- split(x$average.coverage, x$chr)
    
    # for small panels the median appears more robust.
    if (length(xx[[sex.chr[2]]]) < 10) {    
        avg.coverage <- sapply(xx, median, na.rm=TRUE)
    } else {
        if (remove.outliers) xx <- lapply(xx, .removeOutliers)
        avg.coverage <- sapply(xx, mean, na.rm=TRUE)
    }    

    if (is.na(avg.coverage[sex.chr[1]]) || is.na(avg.coverage[sex.chr[2]]) ) {
        flog.warn(
            "Allosome coverage missing, cannot determine sex.")
        return(NA)
    }    
    
    avg.autosome.coverage <- mean(avg.coverage[-match(sex.chr, 
        names(avg.coverage))],na.rm=TRUE)
    autosome.ratio <- avg.autosome.coverage/(avg.coverage[sex.chr[1]]+0.0001)
    if (autosome.ratio > 5) { 
        flog.info("Allosome coverage very low, cannot determine sex.")
        return(NA)
    }
    XY.ratio <- avg.coverage[sex.chr[1]]/ (avg.coverage[sex.chr[2]]+ 0.0001)
    flog.info("Mean coverages: chrX: %.2f, chrY: %.2f, chr1-22: %.2f.",
            avg.coverage[sex.chr[1]], avg.coverage[sex.chr[2]],
            avg.autosome.coverage)
    if (XY.ratio > min.ratio) return("F")
    if (XY.ratio > min.ratio.na) return(NA)
    return("M") 
}

.getSexChr <- function(chrom) {
    if ("chrX" %in% chrom) {
        return(c("chrX", "chrY"))
    } else if ("X" %in% chrom) {
        return(c("X", "Y"))
    }
    return(as.character(23:24))    
}



#' Get sample sex from a VCF file
#' 
#' This function detects non-random distribution of homozygous variants on
#' chromosome X compared to all other chromosomes. A non-significant Fisher's
#' exact p-value indicates more than one chromosome X copy. This function is
#' called in runAbsoluteCN as sanity check when a VCF is provided. It is also
#' useful for determining sex when no sex marker genes on chrY (e.g. AMELY) are
#' available.
#' 
#' 
#' @param vcf CollapsedVCF object, read in with the \code{readVcf} function
#' from the VariantAnnotation package.
#' @param tumor.id.in.vcf The tumor id in the CollapsedVCF (optional).
#' @param min.or Minimum odds-ratio to call sample as male. If p-value is not
#' significant due to a small number of SNPs on chromosome X, sample will be
#' called as NA even when odds-ratio exceeds this cutoff.
#' @param min.or.na Minimum odds-ratio to not call a sample. Odds-ratios in the
#' range \code{min.or.na} to \code{min.or} define a grey area in which samples
#' are not called. Contamination can be a source of ambiguous calls.
#' @param max.pv Maximum Fisher's exact p-value to call sample as male.
#' @param homozygous.cutoff Minimum allelic fraction to call position
#' homozygous.
#' @param af.cutoff Remove all SNVs with allelic fraction lower than the
#' specified value.
#' @param use.somatic.status If somatic status and germline data is available,
#' then exclude somatic variants.
#' @return Returns a \code{character(1)} with \code{M} for male, \code{F} for
#' female, or \code{NA} if unknown.
#' @author Markus Riester
#' @seealso \code{\link{getSexFromCoverage}}
#' @examples
#' 
#' vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
#' vcf <- readVcf(vcf.file, "hg19")
#' # This example vcf is already filtered and contains no homozygous calls,
#' # which are necessary for determining sex from chromosome X.
#' getSexFromVcf(vcf)
#' 
#' @export getSexFromVcf
#' @importFrom stats fisher.test
getSexFromVcf <- function(vcf, tumor.id.in.vcf=NULL, min.or = 4, 
    min.or.na = 2.5, max.pv = 0.001, homozygous.cutoff = 0.95,
    af.cutoff = 0.2, use.somatic.status=TRUE) {
    if (is.null(tumor.id.in.vcf)) {
        tumor.id.in.vcf <- .getTumorIdInVcf(vcf) 
    }
    sex.chr <- .getSexChr(seqlevels(vcf))

    chrY <- seqnames(vcf) == sex.chr[2]
    vcf <- vcf[!chrY]

    if (!is.null(info(vcf)$SOMATIC) && use.somatic.status) {
        vcf <- vcf[!info(vcf)$SOMATIC]
    } else {
        af <- geno(vcf)$FA[,tumor.id.in.vcf] > af.cutoff
        vcf <- vcf[af]
    }

    if (!nrow(vcf)) {
        flog.info("No germline variants in VCF.")
        return(NA)
    }    

    chrX <- seqnames(vcf) == sex.chr[1]

    homozygous <- geno(vcf)$FA[,tumor.id.in.vcf] > homozygous.cutoff
    if ( sum(homozygous)/length(homozygous) < 0.001 ) {
        flog.info("No homozygous variants in VCF, provide unfiltered VCF.")
        return(NA)
    }

    if (!sum(chrX)) {
        flog.info("No variants on chrX in VCF.")
        return(NA)
    }    
    res <- fisher.test(homozygous, as.vector(chrX))
    flog.info("%i homozygous and %i heterozygous variants on chrX.", 
        sum( homozygous & as.vector(chrX)), 
        sum( !homozygous & as.vector(chrX)))
    sex <- "F"    
    if (res$estimate >= min.or.na) sex <- NA
    if (res$estimate >= min.or && res$p.value > max.pv) sex <- NA
    if (res$p.value <= max.pv && res$estimate >= min.or) sex <- "M"
    flog.info("Sex from VCF: %s (Fisher's p-value: %s, odds-ratio: %.2f).", 
        sex, ifelse(res$p.value < 0.0001, "< 0.0001", round(res$p.value, digits=3)), 
        res$estimate)
    return(sex)    
}

.getSex <- function(sex, normal, tumor) {
    if (sex != "?") return(sex)
    sex.tumor <- getSexFromCoverage(tumor)
    sex.normal <- getSexFromCoverage(normal)
    if (!identical(sex.tumor, sex.normal)) {
        flog.warn("Sex tumor/normal mismatch: tumor = ", sex.tumor, 
            " normal = ", sex.normal)
    }
    sex <- sex.tumor    
    if (is.na(sex)) sex = "?"
    sex
}
.fixAllosomeCoverage <- function(sex, tumor) {
    sex.chr <- .getSexChr(tumor$chr)
    if (sex=="M" || sex=="?") {
        tumor <- .removeChr(tumor, remove.chrs=sex.chr)
    } else if (sex=="F") {
        tumor <- .removeChr(tumor, remove.chrs=sex.chr[2])
    }       
    tumor
}

.fixAllosomeSegmentation <- function(sex, seg) {
    if (sex=="diploid") return(seg)
    sex.chr <- .getSexChr(seg$chrom)
    remove.chrs <- sex.chr
    if (sex=="F") {
        remove.chrs <-sex.chr[2]
    }       
    seg[!seg$chrom %in% remove.chrs,]
}

