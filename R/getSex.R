getSexFromCoverage <- structure(function(# Get sample sex from coverage
### This function determines the sex of a sample by the coverage 
### ratio of chrX and chrY. Loss of chromosome Y (LOY) can result in a wrong
### female call. 
gatk.coverage, 
### GATK coverage file or data read with readCoverageGatk.
min.ratio=20,
### Min chrX/chrY coverage ratio to call sample as female.
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
    avg.coverage <- sapply(xx, mean, na.rm=TRUE)
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
    if (XY.ratio > min.ratio) return("F")
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
vcf,
### CollapsedVCF object, read in with the readVcf function 
### from the VariantAnnotation package.
tumor.id.in.vcf=NULL, 
### The tumor id in the CollapsedVCF (optional).
min.or=4,
### Minimum odds-ratio to call sample as male.
min.or.na=2,
### Minimum odds-ratio to not call a sample.
max.pv=0.0001,
### Maximum Fisher's exact p-value to call sample as male.
verbose=TRUE
### Verbose output.
) {
    if (is.null(tumor.id.in.vcf)) {
        tumor.id.in.vcf <- names( which.min(colSums(geno(vcf)$GT=="0")) )
    }
    chrY <- seqnames(vcf) == "chrY" | seqnames(vcf) == "24"
    vcf <- vcf[!chrY]

    chrX <- seqnames(vcf) == "chrX" | seqnames(vcf) == "23"
    homozygous <- geno(vcf)$FA[,tumor.id.in.vcf] == 1
    if ( sum(homozygous)/length(homozygous) < 0.001 ) {
        if (verbose) message("No homozygous variants in VCF, provide unfiltered VCF.")
        return(NA)
    }
    res <- fisher.test(homozygous, as.vector(chrX))
    sex <- "F"    
    if (res$p.value <= max.pv && res$estimate >= min.or.na) sex <- NA
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
