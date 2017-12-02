#' Find best normal sample in database
#' 
#' Function to find the best matching normal for a provided tumor sample.
#' 
#' 
#' @param tumor.coverage.file Coverage file or data of a tumor sample.
#' @param normalDB Database of normal samples, created with
#' \code{\link{createNormalDatabase}}.
#' @param pcs Principal components to use for distance calculation.
#' @param num.normals Return the \code{num.normals} best normals.
#' @param ignore.sex If \code{FALSE}, detects sex of sample and returns best
#' normals with matching sex.
#' @param sex Sex of sample. If \code{NULL}, determine with
#' \code{\link{getSexFromCoverage}} and default parameters. Valid values are
#' \code{F} for female, \code{M} for male. If all chromosomes are diploid,
#' specify \code{diploid}.
#' @param normal.coverage.files Only consider these normal samples. If
#' \code{NULL}, use all in the database. Must match
#' \code{normalDB$normal.coverage.files}.
#' @param pool If \code{TRUE}, use \code{\link{poolCoverage}} to pool best
#' normals.
#' @param pool.weights Either find good pooling weights by optimization or
#' weight all best normals equally.
#' @param plot.pool Allows the pooling function to create plots.
#' @param \dots Additional arguments passed to \code{\link{poolCoverage}}.
#' @return Filename of the best matching normal.
#' @author Markus Riester
#' @seealso \code{\link{createNormalDatabase} \link{getSexFromCoverage}}
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
#' best.normal.coverage.file <- findBestNormal(tumor.coverage.file, normalDB)
#' 
#' pool <- findBestNormal(tumor.coverage.file, normalDB, num.normals=2, 
#'     pool=TRUE)
#' 
#' @export findBestNormal
#' @importFrom stats dist predict
findBestNormal <- function(tumor.coverage.file, normalDB, pcs=1:3, 
    num.normals = 1, ignore.sex = FALSE, sex = NULL, 
    normal.coverage.files = NULL, pool = FALSE, 
    pool.weights = c("voom", "equal"), plot.pool = FALSE, 
    ...) {
    if (is.character(tumor.coverage.file)) {
        tumor  <- readCoverageFile(tumor.coverage.file)
    } else {
        tumor <- tumor.coverage.file
    }    
    if (!is.null(normalDB$gatk.normal.files)) {
        warning("Normal database generated with PureCN version pre 1.2. ",
         "Please regenerate. It will stop working with 1.6.") 
        normalDB$normal.coverage.files <- normalDB$gatk.normal.files
    }    
    normal <- readCoverageFile(normalDB$normal.coverage.files[1])
    if (!identical(as.character(normal) , as.character(tumor))) {
        .stopUserError("tumor.coverage.file and normalDB do not align.")
    }    

    x <- matrix(tumor[normalDB$exons.used]$average.coverage,nrow=1)
    x[is.na(x)] <- 0

    idx.pcs <- pcs
    idx.pcs <- idx.pcs[idx.pcs %in% seq_len(ncol(normalDB$pca$x))]

    idx.normals <- seq_along(normalDB$normal.coverage.files)

    if (!ignore.sex && !is.null(normalDB$sex) && 
        sum(!is.na(normalDB$sex))>0) {
        if (is.null(sex)) {
            sex <- getSexFromCoverage(tumor)
        }
        flog.info("Sample sex: %s", sex)
        if (!is.na(sex)) {
            idx.normals <- which(normalDB$sex == sex)
        }
        if (length(idx.normals) < 2) {
            flog.warn("Not enough samples of sex %s %s", sex, 
                "in database. Ignoring sex.")
            idx.normals <- seq_along(normalDB$normal.coverage.files)
        }
    }

    if (!is.null(normal.coverage.files)) {
        idx.normals <- idx.normals[normalDB$normal.coverage.files[idx.normals] %in%
            normal.coverage.files]
    }    

    best.match <- order(sapply(idx.normals, function(i) 
        dist(rbind(predict(normalDB$pca, x)[1,idx.pcs], 
             predict(normalDB$pca)[i,idx.pcs]))[1]
    ))

    normal.coverage.files <- normalDB$normal.coverage.files[idx.normals][head(best.match, num.normals)]
    if (pool) {
      normals <- lapply(normal.coverage.files, readCoverageFile)
      pool.weights <- match.arg(pool.weights)
      flog.info("Pooling %s.", paste(basename(normal.coverage.files), 
        collapse=", "))
      w <- NULL
      if (pool.weights == "voom" && num.normals > 1) {
          vmlogRatio <- .voomLogRatio(tumor, 
            normal.coverage.files=normal.coverage.files, plot.voom=plot.pool)
          fakeNormal <- tumor
          fakeNormal$average.coverage <- 2 ^ (log2(tumor$average.coverage) - vmlogRatio$logRatio)
          fakeNormal$coverage <- fakeNormal$average.coverage * width(fakeNormal)
          fakeNormal$se <- vmlogRatio$logRatioSe
          return(fakeNormal)
      }
      return(poolCoverage(normals, w=w, ...))
    } 
    normal.coverage.files
}


#' Plot the PCA of tumor and its best normal(s)
#' 
#' This function can be used to understand how a best normal is chosen by the
#' \code{\link{findBestNormal}} function. It can be also used to tune the best
#' normal selection by finding good parameter values for \code{num.normals} and
#' \code{pcs}.
#' 
#' 
#' @param normal.coverage.files Coverage file names of normal samples, 
#' typically identified via \code{\link{findBestNormal}}.
#' @param tumor.coverage.file Coverage file or data  of a tumor sample.
#' @param normalDB Database of normal samples, created with
#' \code{\link{createNormalDatabase}}.
#' @param x Principal component (PC) to be plotted on x-axis.
#' @param y PC to be plotted on y-axis.
#' @param col.tumor Color of tumor in plot.
#' @param col.best.normal Color of best normals in plot.
#' @param col.other.normals Color of other normals in plot.
#' @param \dots Arguments passed to the \code{plot} function.
#' @return Returns \code{NULL}.
#' @author Markus Riester
#' @seealso \code{\link{createNormalDatabase} \link{findBestNormal}}
#' @examples
#' 
#' normal.coverage.file <- system.file("extdata", "example_normal.txt", 
#'     package="PureCN")
#' normal2.coverage.file <- system.file("extdata", "example_normal2.txt",
#'      package="PureCN")
#' normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)
#' normalDB <- createNormalDatabase(normal.coverage.files)
#' 
#' tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
#'     package="PureCN")
#' best.normal.coverage.file <- findBestNormal(tumor.coverage.file, normalDB)
#' plotBestNormal(best.normal.coverage.file, tumor.coverage.file, normalDB)
#' 
#' # Display sample sex. The first point in the plot is always tumor.
#' plotBestNormal(best.normal.coverage.file, tumor.coverage.file, normalDB,  
#'     pch=c(1,ifelse(normalDB$sex=="F", 1, 2)))
#' 
#' @export plotBestNormal
plotBestNormal <- function(normal.coverage.files, tumor.coverage.file,
    normalDB, x = 1, y = 2, col.tumor = "red", col.best.normal = "blue",
    col.other.normals = "black", ...) {
    if (is.character(tumor.coverage.file)) {
        tumor  <- readCoverageFile(tumor.coverage.file)
    } else {
        tumor <- tumor.coverage.file
    }    
    xx <- matrix(tumor[normalDB$exons.used]$average.coverage,nrow=1)
    xx[is.na(xx)] <- 0

    xx <- predict(normalDB$pca,xx)
    xx <- rbind(xx, normalDB$pca$x)
    plot(xx[,x],xx[,y],col=c(col.tumor, 
        ifelse(normalDB$normal.coverage.files %in% 
        normal.coverage.files, col.best.normal, col.other.normals)),
        xlab=paste("PC",x), ylab=paste("PC",y),...)
}
