findBestNormal <- structure(function(#Find best normal sample in database
### Function to find the best matching normal for a provided 
### tumor sample.
tumor.coverage.file, 
### GATK coverage file of a tumor sample.
normalDB,
### Database of normal samples, created with 
### \code{\link{createNormalDatabase}}.
##seealso<< \code{\link{createNormalDatabase} \link{getSexFromCoverage}}
pcs=1:3,
### Principal components to use for distance calculation.
num.normals=1,
### Return the \code{num.normals} best normals.
ignore.sex=FALSE,
### If \code{FALSE}, detects sex of sample and returns best normals
### with matching sex.
sex=NULL,
### Sex of sample. If \code{NULL}, determine with 
### \code{\link{getSexFromCoverage}} and default parameters.
### Valid values are \code{F} for female, \code{M} for male. If all 
### chromosomes are diploid, specify \code{diploid}. 
normal.coverage.files=NULL,
### Only consider these normal samples. If \code{NULL}, use all in 
### the database. Must match \code{normalDB$normal.coverage.files}. 
verbose=TRUE
### Verbose output.
) {
    if (is.character(tumor.coverage.file)) {
        tumor  <- readCoverageGatk(tumor.coverage.file)
    } else {
        tumor <- tumor.coverage.file
    }    
    if (!is.null(normalDB$gatk.normal.files)) {
        warning("Normal database generated with PureCN version pre 1.2. ",
         "Please regenerate. It will stop working with 1.4.") 
        normalDB$normal.coverage.files <- normalDB$gatk.normal.files
    }    
    x <- t(tumor[normalDB$exons.used,"average.coverage", drop=FALSE])
    x[is.na(x)] <- 0

    idx.pcs <- pcs
    idx.pcs <- idx.pcs[idx.pcs %in% seq_len(ncol(normalDB$pca$x))]

    idx.normals <- seq_along(normalDB$normal.coverage.files)

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
            idx.normals <- seq_along(normalDB$normal.coverage.files)
        }
    }

    if (!is.null(normal.coverage.files)) {
        idx.normals <- idx.normals[normalDB$normal.coverage.files[idx.normals] %in%
            normal.coverage.files]
    }    

    best.match <- order(sapply(idx.normals, function(i) 
        dist( rbind(predict(normalDB$pca, x)[1,idx.pcs], 
              predict(normalDB$pca)[i,idx.pcs]))[1]
    ))

    normalDB$normal.coverage.files[idx.normals][head(best.match, num.normals)]
### Filename of the best matching normal.
},ex=function() {
normal.coverage.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
normal2.coverage.file <- system.file("extdata", "example_normal2.txt", 
    package="PureCN")
normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)
normalDB <- createNormalDatabase(normal.coverage.files)

tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")
best.normal.coverage.file <- findBestNormal(tumor.coverage.file, normalDB)
})    


plotBestNormal <- structure(
    function(#Plot the PCA of tumor and its best normal(s)
### This function can be used to understand how a best normal is chosen
### by the \code{\link{findBestNormal}} function. It can be also used 
### to tune the best normal selection by finding good parameter values for
### \code{num.normals} and \code{pcs}.
normal.coverage.files,
### GATK coverage file of normal files, typically identified via 
### \code{\link{findBestNormal}}.
tumor.coverage.file,
### GATK coverage file of a tumor sample.
normalDB,
### Database of normal samples, created with 
### \code{\link{createNormalDatabase}}.
##seealso<< \code{\link{createNormalDatabase} \link{findBestNormal}}
x=1,
### Principal component (PC) to be plotted on x-axis.
y=2,
### PC to be plotted on y-axis.
col.tumor="red",
### Color of tumor in plot.
col.best.normal="blue",
### Color of best normals in plot.
col.other.normals="black",
### Color of other normals in plot.
...
### Arguments passed to the \code{plot} function.
) {
    if (is.character(tumor.coverage.file)) {
        tumor  <- readCoverageGatk(tumor.coverage.file)
    } else {
        tumor <- tumor.coverage.file
    }    
    xx <- t(tumor[normalDB$exons.used,"average.coverage", drop=FALSE])
    xx[is.na(xx)] <- 0

    xx <- predict(normalDB$pca,xx)
    xx <- rbind(xx, normalDB$pca$x)
    plot(xx[,x],xx[,y],col=c(col.tumor, 
        ifelse( normalDB$normal.coverage.files %in% 
        normal.coverage.files, col.best.normal, col.other.normals)),
        xlab=paste("PC",x), ylab=paste("PC",y),...)
### Returns \code{NULL}.
},ex=function() {
normal.coverage.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
normal2.coverage.file <- system.file("extdata", "example_normal2.txt",
     package="PureCN")
normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)
normalDB <- createNormalDatabase(normal.coverage.files)

tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")
best.normal.coverage.file <- findBestNormal(tumor.coverage.file, normalDB)
plotBestNormal(best.normal.coverage.file, tumor.coverage.file, normalDB)

# Display sample sex. The first point in the plot is always tumor.
plotBestNormal(best.normal.coverage.file, tumor.coverage.file, normalDB,  
    pch=c(1,ifelse(normalDB$sex=="F", 1, 2)))
})

