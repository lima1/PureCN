findBestNormal <- structure(function(#Find best normal sample in database
### Function to find the best matching normal for a provided 
### tumor sample.
gatk.tumor.file, 
### GATK coverage file of a tumor sample.
normalDB,
### Database of normal samples, created with 
### createNormalDatabase().
pcs=1:3,
### Principal components to use for distance calculation.
num.normals=1,
### Return the num.normals best normals.
ignore.sex=FALSE,
### If FALSE, detects sex of sample and returns best normals
### with matching sex.
sex=NULL,
### Sex of sample. If NULL, determine with getSexFromCoverage
### and default parameters.
### Valid values are "F" for female, "M" for male. If all 
### chromosomes are diploid, specify "diploid". 
gatk.normal.files=NULL,
### Only consider these normal samples. If NULL, use all in 
### the database. Must match normalDB$gatk.normal.files. 
verbose=TRUE
### Verbose output.
) {
    if (is.character(gatk.tumor.file)) {
        tumor  <- readCoverageGatk(gatk.tumor.file)
    } else {
        tumor <- gatk.tumor.file
    }    
    x <- t(tumor[normalDB$exons.used,"average.coverage", drop=FALSE])
    x[is.na(x)] <- 0

    idx.pcs <- pcs
    idx.pcs <- idx.pcs[idx.pcs %in% seq_len(ncol(normalDB$pca$x))]

    idx.normals <- seq_along(normalDB$gatk.normal.files)

    if (!ignore.sex && !is.null(normalDB$sex) && 
        sum(!is.na(normalDB$sex))>0) {
        if (is.null(sex)) {
            sex <- getSexFromCoverage(tumor, verbose=FALSE)
        }
        if (verbose) message("Sex of sample: ", sex)
        if (!is.na(sex)) {
            idx.normals <- which(normalDB$sex == sex)
        }
        if (length(idx.normals) < 2) { 
            warning("Not enough samples of sex ", sex, 
                " in database. Ignoring sex.")
            idx.normals <- seq_along(normalDB$gatk.normal.files)
        }
    }

    if (!is.null(gatk.normal.files)) {
        idx.normals <- idx.normals[normalDB$gatk.normal.files[idx.normals] %in%
            gatk.normal.files]
    }    

    best.match <- order(sapply(idx.normals, function(i) 
        dist( rbind(predict(normalDB$pca, x)[1,idx.pcs], 
              predict(normalDB$pca)[i,idx.pcs]))[1]
    ))

    normalDB$gatk.normal.files[idx.normals][head(best.match, num.normals)]
### Filename of the best matching normal.
},ex=function() {
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gatk.normal2.file <- system.file("extdata", "example_normal2.txt", 
    package="PureCN")
gatk.normal.files <- c(gatk.normal.file, gatk.normal2.file)
normalDB <- createNormalDatabase(gatk.normal.files)

gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")
gatk.best.normal.file <- findBestNormal(gatk.tumor.file, normalDB)
})    


plotBestNormal <- structure(
    function(#Plot the PCA of tumor and its best normal(s)
### This function can be used to understand how a best normal is chosen
### by the findBestNormal function. It can be also used to tune the
### best normal selection by finding good parameter values for
### num.normals and pcs.
gatk.normal.files,
### GATK coverage file of normal files, typically identified via 
### findBestNormal.
gatk.tumor.file,
### GATK coverage file of a tumor sample.
normalDB,
### Database of normal samples, created with createNormalDatabase().
x=1,
### PC to be plotted on x-axis.
y=2,
### PC to be plotted on y-axis.
col.tumor="red",
### Color of tumor in plot.
col.best.normal="blue",
### Color of best normals in plot.
col.other.normals="black",
### Color of best normals in plot.
...
### Arguments passed to the plot function.
) {
    if (is.character(gatk.tumor.file)) {
        tumor  <- readCoverageGatk(gatk.tumor.file)
    } else {
        tumor <- gatk.tumor.file
    }    
    xx <- t(tumor[normalDB$exons.used,"average.coverage", drop=FALSE])
    xx[is.na(xx)] <- 0

    xx <- predict(normalDB$pca,xx)
    xx <- rbind(xx, normalDB$pca$x)
    plot(xx[,x],xx[,y],col=c(col.tumor, 
        ifelse( normalDB$gatk.normal.files %in% 
        gatk.normal.files, col.best.normal, col.other.normals)),
        xlab=paste("PC",x), ylab=paste("PC",y),...)
### Returns NULL
},ex=function() {
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gatk.normal2.file <- system.file("extdata", "example_normal2.txt",
     package="PureCN")
gatk.normal.files <- c(gatk.normal.file, gatk.normal2.file)
normalDB <- createNormalDatabase(gatk.normal.files)

gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")
gatk.best.normal.file <- findBestNormal(gatk.tumor.file, normalDB)
plotBestNormal(gatk.best.normal.file, gatk.tumor.file, normalDB)

# Display sample sex. The first point in the plot is always tumor.
plotBestNormal(gatk.best.normal.file, gatk.tumor.file, normalDB,  
    pch=c(1,ifelse(normalDB$sex=="F", 1, 2)))
})

