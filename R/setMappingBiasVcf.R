#' Set Mapping Bias VCF
#' 
#' Function to set mapping bias for each variant in the provided
#' \code{CollapsedVCF} object. By default, it returns the same value for all
#' variants, but a pool of normal samples can be provided for position-specific
#' mapping bias calculation.
#' 
#' 
#' @param vcf \code{CollapsedVCF} object, read in with the \code{readVcf}
#' function from the VariantAnnotation package.
#' @param tumor.id.in.vcf Id of tumor in case multiple samples are stored in
#' VCF.
#' @param normal.panel.vcf.file Combined VCF file of a panel of normals,
#' reference and alt counts as AD genotype field. Should be compressed and
#' indexed with bgzip and tabix, respectively. One can provide a precomputed
#' mapping bias database 
#' (obtained by \code{\link{calculateMappingBiasVcf}}).
#' @param min.normals Minimum number of normals with heterozygous SNP for
#' calculating position-specific mapping bias. Requires
#' \code{normal.panel.vcf.file}.
#' @param smooth Impute mapping bias of variants not found in the panel by
#' smoothing of neighboring SNPs. Requires \code{normal.panel.vcf.file}.
#' @param smooth.n Number of neighboring variants used for smoothing.
#' @return A \code{data.frame} with elements \item{bias}{A \code{numeric(nrow(vcf))} 
#' vector with the mapping bias of for each
#' variant in the \code{CollapsedVCF}. Mapping bias is expected as scaling
#' factor. Adjusted allelic fraction is (observed allelic fraction)/(mapping
#' bias). Maximum scaling factor is 1 and means no bias.}
#' \item{pon.count}{A \code{numeric(nrow(vcf))} vector with the number
#' of hits in the \code{normal.panel.vcf.file}.}
#' \item{shape1, shape2}{Fit of a beta distribution.}
#' @author Markus Riester
#' @examples
#' 
#' # This function is typically only called by runAbsoluteCN via 
#' # fun.setMappingBiasVcf and args.setMappingBiasVcf.
#' vcf.file <- system.file("extdata", "example.vcf.gz", package="PureCN")
#' vcf <- readVcf(vcf.file, "hg19")
#' vcf.bias <- setMappingBiasVcf(vcf)        
#' 
#' @export setMappingBiasVcf
setMappingBiasVcf <- function(vcf, tumor.id.in.vcf = NULL,
normal.panel.vcf.file = NULL, min.normals = 2, smooth = TRUE, smooth.n = 5) {

    if (is.null(tumor.id.in.vcf)) {
        tumor.id.in.vcf <- .getTumorIdInVcf(vcf)
    }
    if (min.normals < 2) .stopUserError("min.normals must be >=2.")

    mappingBias <- 1
    if (!is.null(info(vcf)$SOMATIC) && ncol(vcf) > 1) {
         normal.id.in.vcf <- .getNormalIdInVcf(vcf, tumor.id.in.vcf)
         faAll <- as.numeric(geno(vcf)$FA[!info(vcf)$SOMATIC, normal.id.in.vcf])
         mappingBias <- mean(faAll, na.rm=TRUE) * 2
         flog.info("Found SOMATIC annotation in VCF. Setting mapping bias to %.3f.",
            mappingBias)
    }
    if (is.null(info(vcf)$SOMATIC) && is.null(normal.panel.vcf.file)) {
        flog.info(
            "VCF does not contain somatic status. For best results, consider%s%s",
            " providing normal.panel.vcf.file when matched normals are not ",
            "available.")
    }
    tmp <- rep(mappingBias, nrow(vcf))
    # Defines the maximum value for the mapping bias scaling factor.
    # 1 assumes that the reference allele can never have
    # a lower mappability than the alt allele.
    max.bias <- 1.2
    tmp[tmp > max.bias] <- max.bias
    if (is.null(normal.panel.vcf.file)) {
        return(data.frame(bias = tmp, shape1 = NA, shape2 = NA))
    }

    if (file_ext(normal.panel.vcf.file) == "rds") {
        mappingBias <- readRDS(normal.panel.vcf.file)
    } else {
        nvcf <- .readNormalPanelVcfLarge(vcf, normal.panel.vcf.file)
        if (nrow(nvcf) < 1) {
            flog.warn("setMappingBiasVcf: no hits in %s.", normal.panel.vcf.file)
            return(data.frame(bias = tmp, shape1 = NA, shape2 = NA))
        }
        mappingBias <- .findMaxBetaShape(.calculateMappingBias(nvcf, min.normals))
    }
    .annotateMappingBias(tmp, vcf, mappingBias, max.bias, smooth, smooth.n)
}

.annotateMappingBias <- function(tmp, vcf, mappingBias, max.bias, smooth, smooth.n) {
    .extractBias <- function(i, start.var) {
        start <- max(1, i-smooth.n)
        end <- min(length(mappingBias), (i+smooth.n))
        bias <- mappingBias[seq(start,end)]
        #make sure
        bias <- bias[seqnames(bias)==seqnames(mappingBias)[i]]
        #weighted.mean(bias$bias, w = sqrt(bias$pon.count)/sqrt(abs(start(bias) - start.var)+1))
        weighted.mean(bias$bias, w = bias$pon.count)
    }
    ponCnt <- integer(length(tmp))
    shape1 <- rep(NA, length(tmp))
    shape2 <- rep(NA, length(tmp))
    ov <- findOverlaps(vcf, mappingBias, select = "first")
    idx <- !is.na(ov)
    tmp[idx] <- mappingBias$bias[ov][idx]
    ponCnt[idx] <- mappingBias$pon.count[ov][idx]
    if (!is.null(mappingBias$shape1)) {
        shape1[idx] <- mappingBias$shape1[ov][idx]
        shape2[idx] <- mappingBias$shape2[ov][idx]
    }
    if (smooth) {
        flog.info("Imputing mapping bias for %i variants...", 
            sum(!idx, na.rm = TRUE))
        near <- nearest(vcf[!idx], mappingBias, ignore.strand=TRUE)
        start.var <- start(vcf)[!idx]
        tmp[!idx][!is.na(near)] <- sapply(which(!is.na(near)), function(i) .extractBias(near[i], start.var[i]))
    }
    if (anyNA(tmp)) {
        flog.warn("Could not impute mapping bias for all variants. Did you use calculateMappingBiasVcf?")
        tmp[is.na(tmp)] <- 1
    }
    tmp[tmp > max.bias] <- max.bias
    return(data.frame(bias = tmp, pon.count = ponCnt,
                shape1 = shape1, shape2 = shape2))
}
    
.calculateMappingBias <- function(nvcf, min.normals, min.normals.betafit = 7,
                                  min.median.coverage.betafit = 5,
                                  betafit.coverage.outliers = c(0.25, 4)) {
    if (ncol(nvcf) < 2) {
        .stopUserError("The normal.panel.vcf.file contains only a single sample.")
    }
    # TODO: deal with tri-allelic sites
    alt <- apply(geno(nvcf)$AD, c(1,2), function(x) x[[1]][2])
    ref <- apply(geno(nvcf)$AD, c(1,2), function(x) x[[1]][1])
    fa  <- apply(geno(nvcf)$AD, c(1,2), function(x) x[[1]][2]/sum(x[[1]]))
    x <- sapply(seq_len(nrow(nvcf)), function(i) {
        idx <- !is.na(fa[i,]) & fa[i,] > 0.05 & fa[i,] < 0.9
        shapes <- c(NA, NA)
        if (!sum(idx) >= min.normals) return(c(0, 0, 0, 0, shapes))
        dp <- alt[i,] + ref[i,] 
        mdp <- median(dp, na.rm = TRUE)   
        idx2 <- idx & dp >= mdp * betafit.coverage.outliers[1] &
                      dp <= mdp * betafit.coverage.outliers[2] 

        if (sum(idx2) >= min.normals.betafit && mdp >= min.median.coverage.betafit) {
            fit <- try(fitdist(fa[i, idx2], "beta"), silent = TRUE)
            if (class(fit) == "try-error") {
                flog.warn("Could not fit beta dist for %s (%s).",
                    as.character(rowRanges(nvcf[i])),
                    paste0(round(fa[i, idx2], digits = 3), collapse=","))
            } else {    
                shapes <- fit$estimate
            }
        }
        c(sum(ref[i,idx]), sum(alt[i,idx]), sum(idx), mean(fa[i, idx]), shapes)
    })
    # Add an average "normal" SNP (average coverage and allelic fraction > 0.4)
    # as empirical prior
    psMappingBias <- .adjustEmpBayes(x[1:4,]) * 2
    ponCntHits <- apply(geno(nvcf)$AD, 1, function(x)
        sum(!is.na(unlist(x))) / 2)

    tmp <- rowRanges(nvcf)
    mcols(tmp) <- NULL
    tmp$bias <- psMappingBias
    tmp$pon.count <- ponCntHits
    tmp$shape1 <- x[5,]
    tmp$shape2 <- x[6,]
    tmp
}
 
.findMaxBetaShape <- function(x) {
    idx <- !is.na(x$shape1) & x$bias > 0.8
    xx <- x[idx]
    max_shape <- pmax(xx$shape1, xx$shape2)
    # reshape outliers with very high shape parameters
    cutoff <- mean(sapply(split(max_shape, xx$pon.count), quantile, p = 0.95))
    flog.info("Setting max shape parameter to %.2f.", cutoff)
    max_shape <- pmax(x$shape1, x$shape2)
    scale <- cutoff / max_shape
    scale[which(scale > 1)] <- 1
    x$shape1 <- x$shape1 * scale
    x$shape2 <- x$shape2 * scale
    x
}
    
.readNormalPanelVcfLarge <- function(vcf, normal.panel.vcf.file,
    max.file.size=1, geno="AD", expand=FALSE) {
    genome <- genome(vcf)[1]
    if (!file.exists(normal.panel.vcf.file)) {
        .stopUserError("normal.panel.vcf.file ", normal.panel.vcf.file,
            " does not exist.")
    }
    if (file.size(normal.panel.vcf.file) / 1000^3 > max.file.size ||
        nrow(vcf) < 1000) {
        flog.info("Scanning %s...", normal.panel.vcf.file)
        nvcf <- readVcf(TabixFile(normal.panel.vcf.file), genome = genome,
            ScanVcfParam(which = rowRanges(vcf), info = NA, fixed = NA,
            geno = geno))
    } else {
        flog.info("Loading %s...", normal.panel.vcf.file)
        nvcf <- readVcf(normal.panel.vcf.file, genome = genome,
            ScanVcfParam(info = NA, fixed = NA, geno = geno))
        nvcf <- subsetByOverlaps(nvcf, rowRanges(vcf))
    }
    if (expand) nvcf <- expand(nvcf)
    nvcf
}

.adjustEmpBayes <- function(x) {
    # get all SNPs without dramatic bias
    xg <- x[, x[4, ] > 0.4, drop = FALSE]
    if (ncol(xg) < 2) {
        flog.warn("All SNPs in the database have significant mapping bias!%s",
            " Check your database.")
        shape1 <- 0
        shape2 <- 0
    } else {   
        # calculate the average number of ref and alt reads per sample
        shape1 <- sum(xg[1, ]) / sum(xg[3, ])
        shape2 <- sum(xg[2, ]) / sum(xg[3, ])
    }
    # add those as empirical bayes estimate to all SNPs
    x[1, ] <- x[1, ] + shape1
    x[2, ] <- x[2, ] + shape2
    # get the alt allelic fraction for all SNPs
    apply(x, 2, function(y) y[2] / sum(head(y,2)))
}
