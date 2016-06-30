createNormalDatabase <- structure(function(#Create database of normal samples
### Function to create a database of normal samples, used to find 
### a good match for tumor copy number normalization.
gatk.normal.files,
### Vector with file names pointing to GATK coverage files 
### of normal samples. 
sex=NULL,
### Vector of sex."F" for female, "M" for male. If all chromosomes are diploid, specify "diploid". 
### If NULL determine from coverage.
...
### Arguments passed to the prcomp function.
) {
    gatk.normal.files <- normalizePath(gatk.normal.files)
    normals <- lapply(gatk.normal.files, readCoverageGatk)
    normals.m <- do.call(cbind, 
        lapply(normals, function(x) x$average.coverage))
    idx <- complete.cases(normals.m)
    normals.pca <- prcomp(t(normals.m[idx,]), ...)
    sex.determined <- sapply(normals,getSexFromCoverage, verbose=is.null(sex))
    if (is.null(sex)) {
        sex <- sex.determined
    } else {   
        if (length(sex) != length(normals)) {
            stop("Length of gatk.normal.files and sex different")
        } 
        idx.sex <- sex %in% c(NA, "F", "M", "diploid")
        sex[!idx.sex] <- NA
        if(sum(!idx.sex)>0) warning("Unexpected values in sex ignored.")   
        for (i in seq_along(sex.determined)) {
            if (!is.na(sex.determined[i]) && sex[i] != "diploid" &&
                sex.determined[i] != sex[i]) {
                warning("Sex mismatch in ", gatk.normal.files[i], 
                    ". Sex provided is ", sex, ", but could be ", 
                    sex.determined[i])
            }    
        }    
    }
    list(
        gatk.normal.files=gatk.normal.files, 
        pca=normals.pca, 
        exons.used=idx, 
        coverage=apply(normals.m, 2, mean, na.rm=TRUE), 
        exon.median.coverage=apply(normals.m,1,median, na.rm=TRUE),
        exon.log2.sd.coverage=apply(log2(normals.m+1),1,sd, na.rm=TRUE),
        sex=sex
    )
### A normal database that can be used in the findBestNormal function to 
### retrieve good matching normal samples for a given tumor sample.
},ex=function() {
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gatk.normal2.file <- system.file("extdata", "example_normal2.txt", 
    package="PureCN")
gatk.normal.files <- c(gatk.normal.file, gatk.normal2.file)
normalDB <- createNormalDatabase(gatk.normal.files)
})    

createExonWeightFile <- structure(function(# Calculate exon weights
### Creates an exon weight file useful for segmentation. Requires a
### set of GATK coverage files from normal samples. A small number of 
### tumor (or other normal) samples is then tested against all normals. Exon
### weights will be set proportional to the inverse of coverage standard
### deviation across all normals. Exons with high variance in coverage in the
### pool of normals are thus down-weighted.
gatk.tumor.files, 
### A small number (1-3) of GATK tumor or normal coverage samples.
gatk.normal.files,
### A large number of GATK normal coverage samples (>20) 
### to estimate exon log-ratio standard deviations.
### Should not overlap with files in gatk.tumor.files.
exon.weight.file
### Output filename.
) {
    tumor.coverage <- lapply(gatk.tumor.files,  readCoverageGatk)
    lrs <- lapply(tumor.coverage, function(tc) sapply(gatk.normal.files, 
            function(x) .calcLogRatio(readCoverageGatk(x), tc, verbose=TRUE)))

    lrs <- do.call(cbind, lrs)

    lrs[is.infinite(lrs)] <- NA

    lrs.sd <- apply(lrs, 1, sd, na.rm=TRUE)
    lrs.cnt.na <- apply(lrs,1, function(x) sum(is.na(x)))

    zz <- quantile(lrs.sd,probs=0.75, na.rm=TRUE)/lrs.sd
    zz[zz>1] <-1
    idx <- is.na(zz) | lrs.cnt.na > ncol(lrs)/3
    zz[idx] <- min(zz, na.rm=TRUE)
    ret <- data.frame(Target=tumor.coverage[[1]][,1], Weights=zz)

    write.table(ret, file=exon.weight.file,row.names=FALSE, quote=FALSE, 
        sep="\t")
    invisible(ret)
###A data.frame with exon weights.
}, ex=function() {
exon.weight.file <- "exon_weights.txt"
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gatk.normal2.file <- system.file("extdata", "example_normal2.txt", 
    package="PureCN")
gatk.normal.files <- c(gatk.normal.file, gatk.normal2.file)
gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")

createExonWeightFile(gatk.tumor.file, gatk.normal.files, exon.weight.file)
})
