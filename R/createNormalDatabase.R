#' Create database of normal samples
#' 
#' Function to create a database of normal samples, used to find a good match
#' for tumor copy number normalization. Internally, this function determines
#' the sex of the samples and trains a PCA that is later used for clustering a
#' tumor file with all normal samples in the database.
#' 
#' 
#' @param normal.coverage.files Vector with file names pointing to GATK
#' coverage files of normal samples.
#' @param sex \code{character(length(normal.coverage.files))} with sex for all
#' files.  \code{F} for female, \code{M} for male. If all chromosomes are
#' diploid, specify \code{diploid}. If \code{NULL}, determine from coverage.
#' @param max.mean.coverage Assume that coverages above this value do not
#' necessarily improve copy number normalization. Internally, samples with
#' coverage higher than this value will be normalized to have mean coverage
#' equal to this value. If \code{NULL}, use the 80 percentile as cutoff.
#' If \code{NA}, does not use a maxmimum value.
#' @param \dots Arguments passed to the \code{prcomp} function.
#' @return A normal database that can be used in the
#' \code{\link{findBestNormal}} function to retrieve good matching normal
#' samples for a given tumor sample.
#' @author Markus Riester
#' @seealso \code{\link{findBestNormal}}
#' @examples
#' 
#' normal.coverage.file <- system.file("extdata", "example_normal.txt", 
#'     package="PureCN")
#' normal2.coverage.file <- system.file("extdata", "example_normal2.txt", 
#'     package="PureCN")
#' normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)
#' normalDB <- createNormalDatabase(normal.coverage.files)
#' 
#' @export createNormalDatabase
#' @importFrom stats prcomp
createNormalDatabase <- function(normal.coverage.files, sex = NULL,
max.mean.coverage = NULL, ... ) {
    normal.coverage.files <- normalizePath(normal.coverage.files)
    normals <- .readNormals(normal.coverage.files)

    normals.m <- do.call(cbind, 
        lapply(normals, function(x) x$average.coverage))
    idx <- complete.cases(normals.m)

    z <- apply(normals.m[idx,],2,mean)
    if (is.null(max.mean.coverage)) max.mean.coverage <- 
        quantile(round(z), p=0.8)
    if (!is.na(max.mean.coverage)) { 
        flog.info("Setting maximum coverage in normalDB to %.0f", 
            max.mean.coverage)
        z <- sapply(max.mean.coverage/z, min,1)
        normals.m <- scale(normals.m, 1/z, center=FALSE)
    }
    normals.pca <- prcomp(t(normals.m[idx,]), ...)
    sex.determined <- sapply(normals,getSexFromCoverage)
    if (is.null(sex)) {
        sex <- sex.determined
    } else {   
        if (length(sex) != length(normals)) {
            .stopUserError("Length of normal.coverage.files and sex different")
        } 
        idx.sex <- sex %in% c(NA, "F", "M", "diploid")
        sex[!idx.sex] <- NA
        if(sum(!idx.sex)>0) warning("Unexpected values in sex ignored.")   
        for (i in seq_along(sex.determined)) {
            if (!is.na(sex.determined[i]) && sex[i] != "diploid" &&
                sex.determined[i] != sex[i]) {
                flog.warn("Sex mismatch in %s. Sex provided is %s, but could be %s.", 
                    normal.coverage.files[i], sex[i], sex.determined[i])
            }    
        }    
    }
    list(
        normal.coverage.files=normal.coverage.files, 
        pca=normals.pca, 
        exons.used=idx, 
        coverage=apply(normals.m, 2, mean, na.rm=TRUE), 
        exon.median.coverage=apply(normals.m,1,median, na.rm=TRUE),
        exon.log2.sd.coverage=apply(log2(normals.m+1),1,sd, na.rm=TRUE),
        fraction.missing=apply(normals.m,1,function(x) 
            sum(is.na(x)|x<0.01))/ncol(normals.m),
        sex=sex
    )
}


#' Calculate exon weights
#' 
#' This function is defunct. Please use \code{\link{createTargetWeights}}
#' instead.
#' 
#' 
#' @author Markus Riester
#' @export createExonWeightFile
createExonWeightFile <- function() {
    .Defunct("createTargetWeights")
}



#' Calculate target weights
#' 
#' Creates a target weight file useful for segmentation. Requires a set of GATK
#' coverage files from normal samples. A small number of tumor (or other
#' normal) samples is then tested against all normals. Target weights will be
#' set proportional to the inverse of coverage standard deviation across all
#' normals. Targets with high variance in coverage in the pool of normals are
#' thus down-weighted.
#' 
#' 
#' @param tumor.coverage.files A small number (1-3) of GATK tumor or normal
#' coverage samples.
#' @param normal.coverage.files A large number of GATK normal coverage samples
#' (>20) to estimate target log-ratio standard deviations. Should not overlap
#' with files in \code{tumor.coverage.files}.
#' @param target.weight.file Output filename.
#' @return A \code{data.frame} with target weights.
#' @author Markus Riester
#' @examples
#' 
#' target.weight.file <- "target_weights.txt"
#' normal.coverage.file <- system.file("extdata", "example_normal.txt", 
#'     package="PureCN")
#' normal2.coverage.file <- system.file("extdata", "example_normal2.txt", 
#'     package="PureCN")
#' normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)
#' tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
#'     package="PureCN")
#' 
#' createTargetWeights(tumor.coverage.file, normal.coverage.files, target.weight.file)
#' 
#' @export createTargetWeights
createTargetWeights <- function(tumor.coverage.files, normal.coverage.files,
target.weight.file) {
    flog.info("Loading coverage data...")
    tumor.coverage <- lapply(tumor.coverage.files,  readCoverageGatk)
    normal.coverage <- lapply(normal.coverage.files,  readCoverageGatk)
    lrs <- lapply(tumor.coverage, function(tc) sapply(normal.coverage, 
            function(nc) calculateLogRatio(nc, tc)))

    lrs <- do.call(cbind, lrs)

    lrs[is.infinite(lrs)] <- NA

    lrs.sd <- apply(lrs, 1, sd, na.rm=TRUE)
    lrs.cnt.na <- apply(lrs,1, function(x) sum(is.na(x)))

    zz <- quantile(lrs.sd,probs=0.75, na.rm=TRUE)/lrs.sd
    zz[zz>1] <-1
    idx <- is.na(zz) | lrs.cnt.na > ncol(lrs)/3
    zz[idx] <- min(zz, na.rm=TRUE)
    ret <- data.frame(Target=tumor.coverage[[1]][,1], Weights=zz)

    write.table(ret, file=target.weight.file,row.names=FALSE, quote=FALSE, 
        sep="\t")
    invisible(ret)
}

.readNormals <- function(normal.coverage.files) {
    normals <- lapply(normal.coverage.files, readCoverageGatk)

    # check that all files used the same interval file.
    for (i in seq_along(normals)) {
        if (!identical(normals[[i]]$probe, normals[[1]]$probe)) {
            .stopUserError("All normal.coverage.files must have the same ",
                "intervals. ", normal.coverage.files[i], " is different.")
        }
    }
    normals
}
