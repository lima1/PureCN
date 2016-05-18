getSexFromCoverage <- structure(function(# Get sample sex from coverage
### This function determines the sex of a sample by the coverage 
### ratio of chrX and chrY. Loss of chromosome Y (LOY) can result in a wrong
### female call. 
gatk.coverage, 
### GATK coverage file or data read with readCoverageGatk.
min.ratio=20,
### Min chrX/chrY coverage ratio to call sample as female.
min.ratio.na=7,
### Min chrX/chrY coverage ratio to call sample as NA. This ratio defines a 
### grey zone from min.ratio.na to min.ratio in which samples are not called.
### The default is set to a copy number ratio that would be rare in male samples,
### but lower than expected in female samples. Contamination can be a 
### source of ambiguous calls.
remove.outliers=TRUE,
### Removes coverage outliers before calculating mean chromosome coverages.
verbose=TRUE
### Verbose output.
) {
    if (is.character(gatk.coverage)) {
        x <- readCoverageGatk(gatk.coverage)
    } else {
        x <- gatk.coverage
    }

    sex.chr <- .getSexChr(x)
    xx <- split(x$average.coverage, x$chr)
    
    # for small panels the median appears more robust.
    if (length(xx[[sex.chr[2]]]) < 10) {    
        avg.coverage <- sapply(xx, median, na.rm=TRUE)
    } else {
        if (remove.outliers) xx <- lapply(xx, .removeOutliers)
        avg.coverage <- sapply(xx, mean, na.rm=TRUE)
    }    

    if (is.na(avg.coverage[sex.chr[1]]) || is.na(avg.coverage[sex.chr[2]]) ) {
        if (verbose) message(
            "Allosome coverage appears to be missing, cannot determine sex.")
        return(NA)
    }    

    autosome.ratio <- mean(avg.coverage[-match(sex.chr, names(avg.coverage))], 
        na.rm=TRUE)/(avg.coverage[sex.chr[1]]+0.0001)
    if (autosome.ratio > 5) { 
        if (verbose) message(
            "Allosome coverage very low, cannot determine sex.")
        return(NA)
    }
    XY.ratio <- avg.coverage[sex.chr[1]]/ (avg.coverage[sex.chr[2]]+ 0.0001)
    if (verbose) {
        message("Mean coverage chrX: ",  avg.coverage[sex.chr[1]], 
                ".\nMean coverage chrY: ", avg.coverage[sex.chr[2]])
    }     
    if (XY.ratio > min.ratio) return("F")
    if (XY.ratio > min.ratio.na) return(NA)
    return("M") 
### Returns "M" for male, "F" for female, or NA if unknown.    
}, ex=function(){
    gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
        package="PureCN")
    sex <- getSexFromCoverage(gatk.tumor.file)
})

.getSexChr <- function(gatk.coverage) {
    if ("chrX" %in% gatk.coverage$chr) {
        return(c("chrX", "chrY"))
    }
    return(as.character(23:24))    
}

getSexFromVcf <- structure(function(# Get sample sex from VCF file
### This function detects non-random distribution of homozygous
### variants on chromosome X compared to all other chromosomes.
### A non-significant Fisher's exact p-value indicates a diploid 
### chromosome X.
vcf,
### CollapsedVCF object, read in with the readVcf function 
### from the VariantAnnotation package.
tumor.id.in.vcf=NULL, 
### The tumor id in the CollapsedVCF (optional).
min.or=4,
### Minimum odds-ratio to call sample as male. If p-value is
### not significant due to a small number of SNPs on chromosome X,
### sample will be called as NA.
min.or.na=2.5,
### Minimum odds-ratio to not call a sample. Odds-ratios in the
### range min.or.na to min.or define a grey area in which samples
### are not called. Contamination can be a source of ambiguous calls.
max.pv=0.001,
### Maximum Fisher's exact p-value to call sample as male.
homozygous.cutoff=0.95,
### Minimum allelic fraction to call position homozygous.
af.cutoff=0.03,
### Remove all SNVs with allelic fraction lower than the
### specified value.
verbose=TRUE
### Verbose output.
) {
    if (is.null(tumor.id.in.vcf)) {
        tumor.id.in.vcf <- names( which.min(colSums(geno(vcf)$GT=="0")) )
    }
    chrY <- seqnames(vcf) == "chrY" | seqnames(vcf) == "24"
    vcf <- vcf[!chrY]
    af <- geno(vcf)$FA[,tumor.id.in.vcf] > af.cutoff
    vcf <- vcf[af]

    chrX <- seqnames(vcf) == "chrX" | seqnames(vcf) == "23"
    homozygous <- geno(vcf)$FA[,tumor.id.in.vcf] > homozygous.cutoff
    if ( sum(homozygous)/length(homozygous) < 0.001 ) {
        if (verbose) message("No homozygous variants in VCF, provide unfiltered VCF.")
        return(NA)
    }
    res <- fisher.test(homozygous, as.vector(chrX))
    if (verbose) message(sum( homozygous & as.vector(chrX)), 
        " homozygous and ", sum( !homozygous & as.vector(chrX)), 
        " heterozygous variants on chromosome X.")
    sex <- "F"    
    if (res$estimate >= min.or.na) sex <- NA
    if (res$estimate >= min.or && res$p.value > max.pv) sex <- NA
    if (res$p.value <= max.pv && res$estimate >= min.or) sex <- "M"
    if (verbose) message("Sex from VCF: ", sex, " (Fisher's p-value: ", res$p.value, "  odds-ratio: ", res$estimate, ")")
    return(sex)    
### Returns "M" for male, "F" for female, or NA if unknown.    
}, ex=function() {
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
vcf <- readVcf(vcf.file, "hg19")
# This example vcf is already filtered and contains no homozygous calls,
# which are necessary for determining sex from chromosome X.
getSexFromVcf(vcf)
})    
