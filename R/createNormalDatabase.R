#' Create database of normal samples
#' 
#' Function to create a database of normal samples, used to normalize
#' tumor coverages.
#' 
#' 
#' @param normal.coverage.files Vector with file names pointing to 
#' coverage files of normal samples.
#' @param sex \code{character(length(normal.coverage.files))} with sex for all
#' files.  \code{F} for female, \code{M} for male. If all chromosomes are
#' diploid, specify \code{diploid}. If \code{NULL}, determine from coverage.
#' @param coverage.outliers Exclude samples with coverages below or above
#' the specified cutoffs (fractions of the normal sample coverages median).
#' Only for databases with more than 5 samples.
#' @param min.coverage Exclude intervals with coverage lower than 
#' the specified fraction of the chromosome median in the pool of normals.
#' @param max.missing Exclude intervals with zero coverage in the
#' specified fraction of normal samples.
#' @param low.coverage Specifies the maximum number of total reads 
#' (NOT average coverage) to call a target low coverage.
#' @param plot Diagnostics plot, useful to tune parameters.
#' @param \dots Arguments passed to the \code{prcomp} function.
#' @return A normal database that can be used in the
#' \code{\link{calculateTangentNormal}} function to retrieve a coverage
#' normalization sample for a given tumor sample.
#' @author Markus Riester
#' @seealso \code{\link{calculateTangentNormal}}
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
#' @importFrom Matrix tcrossprod
createNormalDatabase <- function(normal.coverage.files, sex = NULL,
coverage.outliers = c(0.25, 4), min.coverage = 0.25,
max.missing = 0.03, low.coverage = 15, plot = FALSE, ...) {
    normal.coverage.files <- normalizePath(normal.coverage.files)

    if (any(duplicated(normal.coverage.files))) {
        flog.warn("Some normal.coverage.files duplicated.")
        normal.coverage.files <- unique(normal.coverage.files)
    }
    if (length(normal.coverage.files) < 2) {
        .stopUserError("At least 2 normal.coverage.files required.")
    }

    normals <- .readNormals(normal.coverage.files)

    normals.m <- do.call(cbind, 
        lapply(normals, function(x) x$counts))

    low.coverage.targets <- .warnLowCoverageTargets(normals.m, normals[[1]], 
        low.coverage)
    normals.m[is.na(normals.m)] <- 0

    z <- apply(normals.m,2,mean)
    idx.failed <- rep(FALSE, length(normals))

    if (length(normals) > 5) {
        idx.failed <- z < median(z) * coverage.outliers[1] | 
                      z > median(z) * coverage.outliers[2]
        if (sum(idx.failed)) {              
            flog.info("Dropping %s due to outlier coverage.", 
                paste(basename(normal.coverage.files[idx.failed]),collapse=", "))
        }
        normals <- normals[!idx.failed]
        normals.m <- normals.m[,!idx.failed,drop = FALSE]
        normal.coverage.files <- normal.coverage.files[!idx.failed]
    }
    sex.determined <- sapply(normals,getSexFromCoverage)
    if (is.null(sex)) {
        sex <- sex.determined
    } else {
        if (length(sex) != length(idx.failed)) {
            .stopUserError("Length of normal.coverage.files and sex different")
        }
        sex <- sex[!idx.failed]
        idx.sex <- sex %in% c(NA, "F", "M", "diploid")
        sex[!idx.sex] <- NA
        if (sum(!idx.sex)>0) warning("Unexpected values in sex ignored.")
        for (i in seq_along(sex.determined)) {
            if (!is.na(sex.determined[i]) && sex[i] != "diploid" &&
                sex.determined[i] != sex[i]) {
                flog.warn("Sex mismatch in %s. Sex provided is %s, but could be %s.", 
                    normal.coverage.files[i], sex[i], sex.determined[i])
            }    
        }    
    }


    groups <- lapply(c(TRUE, FALSE), function(on.target) {
        idx <- normals[[1]]$on.target == on.target
        intervals <- normals[[1]][idx]
        if (length(intervals)) {
            flog.info("Processing %s-target regions...", 
                ifelse(on.target, "on", "off") )
        }    
        .standardizeNormals(normals.m[idx,], normals[[1]][idx], min.coverage, 
            max.missing, sex)
    })
    
    # merge back some on-target and off-target interval statistics 
    intervals.used <- logical(length(normals[[1]]))
    fraction.missing <- double(length(normals[[1]]))
    intervals.used[normals[[1]]$on.target] <- groups[[1]]$intervals.used
    fraction.missing[normals[[1]]$on.target] <- groups[[1]]$fraction.missing

    if (!is.null(groups[[2]]$intervals.used)) {
        intervals.used[!normals[[1]]$on.target] <- groups[[2]]$intervals.used
        fraction.missing[!normals[[1]]$on.target] <- groups[[2]]$fraction.missing
    }
        
    normalDB <- list(
        normal.coverage.files = normal.coverage.files,
        intervals = as.character(normals[[1]]),
        groups = groups,
        intervals.used = intervals.used,
        sex = sex,
        low.coverage.targets = low.coverage.targets,
        version = 8
    )
    normalDB <- .calculateIntervalWeights(normalDB, normals, plot = plot)
    normalDB
}

.warnLowCoverageTargets <- function(counts, intervals, low.coverage) {
    all.low <- apply(counts, 1, max) <= low.coverage
    low.intervals <- intervals[which(all.low & intervals$on.target)]
    mcols(low.intervals) <- NULL
    n.low <- length(low.intervals)
    if (n.low > 0) {
        flog.info("%i on-target bins with low coverage in all samples.", n.low)
    }
    if (n.low > sum(intervals$on.target, na.rm = TRUE) * 0.05) {
        flog.warn("You are likely not using the correct baits file!")
    }
    low.intervals    
}
     
.standardizeNormals <- function(counts, intervals, min.coverage, max.missing, sex) {
    if (!length(intervals)) return(list(present = FALSE))
    # recalculate without dropped samples
    fcnts <- apply(counts, 2, function(x) x/sum(x))
    fcnts_interval_medians <- apply(fcnts, 1, median, na.rm=TRUE)
    fraction.missing <- apply(counts, 1, function(x)
                    sum(is.na(x)|x<=0))/ncol(counts)

    intervals.used <- .filterIntervalsCreateNormalDB(intervals, 
        fcnts_interval_medians, fraction.missing, min.coverage, max.missing) 

    fcnts <- apply(counts[intervals.used,], 2, function(x) x/sum(x))
    fcnts_interval_medians <- apply(fcnts, 1, median)
    fcnts_interval_medians_F <- fcnts_interval_medians
    fcnts_interval_medians_M <- fcnts_interval_medians
    
    # do we need a sex-aware database?
    if ("F" %in% sex && "M" %in% sex) {
        fcnts_interval_medians_F <- apply(fcnts[, which(sex=="F"), drop = FALSE], 1, median)
        fcnts_interval_medians_M <- apply(fcnts[, which(sex=="M"), drop = FALSE], 1, median)
        sex.chr <- .getSexChr(seqlevels(intervals))
        idx <- as.character(seqnames(intervals)[intervals.used]) %in% sex.chr
        fcnts_interval_medians_F[!idx] <- fcnts_interval_medians[!idx]
        fcnts_interval_medians_M[!idx] <- fcnts_interval_medians[!idx]
        fcnts_std <- fcnts
        for (i in seq(ncol(fcnts))) {
            if (is.null(sex[i]) || is.na(sex[i])) sex[i] <- "?"
            iv <- switch(sex[i],
                "F"=fcnts_interval_medians_F,
                "M"=fcnts_interval_medians_M,
                fcnts_interval_medians)

            fcnts_std[, i] <- fcnts_std[,i] / iv
        }
    } else {  
        fcnts_std <- apply(fcnts,2,function(x) x/fcnts_interval_medians)
    }
    fcnts_interval_non_zero_medians <- apply(fcnts_std, 1, function(x) median(x[x>0]))
    fcnts_std_imp <- apply(fcnts_std, 2, function(x) { x[x<=0] <- fcnts_interval_non_zero_medians[x<=0]; x})
    p=0.001
    li <- quantile(as.vector(fcnts_std_imp),probs=c(p, 1-p))
    fcnts_std_trunc <- fcnts_std_imp
    fcnts_std_trunc[fcnts_std_imp<li[1]] <- li[1]
    fcnts_std_trunc[fcnts_std_imp>li[2]] <- li[2]
    fcnts_std_final <- apply(fcnts_std_trunc, 2, function(x) log2(x/median(x)))
    fcnts_std_final - median(apply(fcnts_std_final,2,median))
    s <- svd(fcnts_std_final)

    list(
        projection = s$u,
        projection.v = s$v,
        intervals.used = intervals.used,
        interval.median.coverage = list(all = fcnts_interval_medians,
                                        F = fcnts_interval_medians_F,        
                                        M = fcnts_interval_medians_M),
        fraction.missing = fraction.missing,
        present=TRUE
    )
}

.denoiseSample <- function(x, normalDB, num.eigen, sex) {
    fcnts <- x$counts[normalDB$intervals.used]
    fcnts <- fcnts/sum(fcnts, na.rm=TRUE)
    if (is.null(sex) || is.na(sex)) sex <- "?"
    iv <- switch(sex,
        "F"=normalDB$interval.median.coverage$F,
        "M"=normalDB$interval.median.coverage$M,
        normalDB$interval.median.coverage$all)

    fcnts_std <- fcnts/iv
    # impute intervals which have no data in sample but in database
    # otherwise the projection will return only NA
    idxMissing <- which(fcnts_std <= 0 | is.na(fcnts_std))
    fcnts_std[idxMissing] <- 1e-9

    fcnts_std_final <- log2(fcnts_std/median(fcnts_std, na.rm = TRUE))
    x$log.ratio.std <- 0.
    x$log.ratio.std[!normalDB$intervals.used] <- NA
    x$log.ratio <- x$log.ratio.std
    x$log.ratio.std[normalDB$intervals.used] <- fcnts_std_final
    if (num.eigen > ncol(normalDB$projection)) num.eigen <- ncol(normalDB$projection)
    P <- normalDB$projection[,seq(1,num.eigen)]    
    x$log.ratio[normalDB$intervals.used] <- fcnts_std_final - 
        as.vector(tcrossprod(fcnts_std_final %*% P, P))
    x
}

#' Calculate tangent normal
#' 
#' Reimplementation of GATK4 denoising. Please cite the relevant GATK
#' publication if you use this in a publication.
#' 
#' 
#' @param tumor.coverage.file Coverage file or data  of a tumor sample.
#' @param normalDB Database of normal samples, created with
#' \code{\link{createNormalDatabase}}.
#' @param num.eigen Number of eigen vectors used.
#' @param ignore.sex If \code{FALSE}, detects sex of sample and returns best
#' normals with matching sex.
#' @param sex Sex of sample. If \code{NULL}, determine with
#' \code{\link{getSexFromCoverage}} and default parameters. Valid values are
#' \code{F} for female, \code{M} for male. If all chromosomes are diploid,
#' specify \code{diploid}.
#' @seealso \code{\link{createNormalDatabase}}
#' @author Markus Riester
#' @examples
#'
#' tumor.coverage.file <- system.file('extdata', 'example_tumor.txt', 
#'     package='PureCN')
#' normal.coverage.file <- system.file("extdata", "example_normal.txt", 
#'     package="PureCN")
#' normal2.coverage.file <- system.file("extdata", "example_normal2.txt", 
#'     package="PureCN")
#' normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)
#' normalDB <- createNormalDatabase(normal.coverage.files)
#' pool <- calculateTangentNormal(tumor.coverage.file, normalDB)
#' 
#' @export calculateTangentNormal
calculateTangentNormal <- function(tumor.coverage.file, normalDB, 
                                   num.eigen = 20, ignore.sex = FALSE, 
                                   sex = NULL) {

    if (is.character(tumor.coverage.file)) {
        tumor  <- readCoverageFile(tumor.coverage.file)
    } else {
        tumor <- tumor.coverage.file
    }
    
    if (!.checkNormalDB(tumor, normalDB)) {
        .stopUserError("tumor.coverage.file and normalDB do not align.")
    }

    if (!ignore.sex && !is.null(normalDB$sex) && 
        sum(!is.na(normalDB$sex))>0) {
        if (is.null(sex)) {
            sex <- getSexFromCoverage(tumor)
        }
        flog.info("Sample sex: %s", sex)
    }
    tumor$log.ratio <- 0.
    tumor$log.ratio.std <- 0.
    denoised <- .denoiseSample(tumor[tumor$on.target], normalDB$groups[[1]], num.eigen, sex)
    ov <- findOverlaps(tumor, denoised)
    tumor$log.ratio[queryHits(ov)] <- denoised$log.ratio[subjectHits(ov)]
    tumor$log.ratio.std[queryHits(ov)] <- denoised$log.ratio.std[subjectHits(ov)]

    if (normalDB$groups[[2]]$present) {
        denoised <- .denoiseSample(tumor[!tumor$on.target], normalDB$groups[[2]], num.eigen, sex)
        ov <- findOverlaps(tumor, denoised)
        tumor$log.ratio[queryHits(ov)] <- denoised$log.ratio[subjectHits(ov)]
        tumor$log.ratio.std[queryHits(ov)] <- denoised$log.ratio.std[subjectHits(ov)]
    }
        
    fakeNormal <- tumor
    idxNA <- is.na(fakeNormal$log.ratio)
    fakeNormal$average.coverage[!idxNA] <- (2 ^ (log2(tumor$average.coverage) - tumor$log.ratio))[!idxNA]
    fakeNormal$coverage <- fakeNormal$average.coverage * width(fakeNormal)
    fakeNormal
}
    
.readNormals <- function(normal.coverage.files) {
    normals <- lapply(normal.coverage.files, readCoverageFile)

    # check that all files used the same interval file.
    for (i in seq_along(normals)) {
        if (!identical(as.character(normals[[i]]), as.character(normals[[1]]))) {
            .stopUserError("All normal.coverage.files must have the same ",
                "intervals. ", normal.coverage.files[i], " is different.")
        }
    }
    normals
}

.filterIntervalsCreateNormalDB <- function(intervals, 
interval.median.coverage, fraction.missing,
min.coverage, max.missing) {
    
    nBefore <- length(intervals)
    intervals.used <- is.finite(fraction.missing) 

    key <- paste(as.character(seqnames(intervals)), intervals$on.target)
    min.coverage <- (sapply(split(interval.median.coverage, 
            key), median, na.rm=TRUE) * min.coverage)[key]

    intervals.used <- intervals.used & !is.na(interval.median.coverage) & 
        interval.median.coverage >= min.coverage
        
    nAfter <- sum(intervals.used)

    if (nAfter < nBefore) {
        flog.info("Removing %i intervals with low coverage in normalDB.", 
            nBefore-nAfter)
    }

    nBefore <- sum(intervals.used)
    intervals.used <- intervals.used & !is.na(fraction.missing) &
        fraction.missing <= max.missing
    nAfter <- sum(intervals.used)

    if (nAfter < nBefore) {
        flog.info("Removing %i intervals with zero coverage in more than %.0f%% of normalDB.", 
            nBefore-nAfter, max.missing*100)
    }

    if (nAfter < 2) {
        .stopUserError("No intervals left after coverage checks.")
    }

    intervals.used
}
