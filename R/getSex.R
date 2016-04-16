getSexFromCoverage <- structure(function(# Get sample sex from coverage
### This function determines the sex of a sample by the coverage 
### ratio of chrX and chrY.
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
