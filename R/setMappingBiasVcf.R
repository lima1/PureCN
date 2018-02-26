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
#' @return A list with elements \item{bias}{A \code{numeric(nrow(vcf))} 
#' vector with the mapping bias of for each
#' variant in the \code{CollapsedVCF}. Mapping bias is expected as scaling
#' factor. Adjusted allelic fraction is (observed allelic fraction)/(mapping
#' bias). Maximum scaling factor is 1 and means no bias.}
#' \item{pon.count}{A \code{numeric(nrow(vcf))} vector with the number
#' of hits in the \code{normal.panel.vcf.file}.}
#' @author Markus Riester
#' @examples
#' 
#' # This function is typically only called by runAbsoluteCN via 
#' # fun.setMappingBiasVcf and args.setMappingBiasVcf.
#' vcf.file <- system.file("extdata", "example_vcf.vcf.gz", package="PureCN")
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
        return(list(bias = tmp))
    }

    if (file_ext(normal.panel.vcf.file) == "rds") {
        mappingBias <- readRDS(normal.panel.vcf.file)
    } else {
        nvcf <- .readNormalPanelVcfLarge(vcf, normal.panel.vcf.file)
        if (nrow(nvcf) < 1) {
            flog.warn("setMappingBiasVcf: no hits in %s.", normal.panel.vcf.file)
            return(list(bias = tmp))
        }
        mappingBias <- .calculateMappingBias(nvcf, min.normals)
    }
    ov <- findOverlaps(vcf, mappingBias, select = "first")

    ponCnt <- integer(length(tmp))

    tmp[!is.na(ov)] <- mappingBias$bias[ov][!is.na(ov)]
    tmp[tmp > max.bias] <- max.bias

    ponCnt[!is.na(ov)] <- mappingBias$pon.count[ov][!is.na(ov)]
    if (smooth) {
        tmpSmoothed <- .smoothVectorByChromosome(tmp,
            as.character(seqnames(vcf)), smooth.n)
        # only smooth variants not found in database
        tmp[is.na(ov)] <- tmpSmoothed[is.na(ov)]
    }
    return(list(bias = tmp, pon.count = ponCnt))
}

.calculateMappingBias <- function(nvcf, min.normals) {
    if (ncol(nvcf) < 2) {
        .stopUserError("The normal.panel.vcf.file contains only a single sample.")
    }
    # TODO: deal with tri-allelic sites
    alt <- apply(geno(nvcf)$AD, c(1,2), function(x) x[[1]][2])
    ref <- apply(geno(nvcf)$AD, c(1,2), function(x) x[[1]][1])
    fa  <- apply(geno(nvcf)$AD, c(1,2), function(x) x[[1]][2]/sum(x[[1]]))
    psMappingBias <- sapply(seq_len(nrow(nvcf)), function(i) {
        idx <- !is.na(fa[i,]) & fa[i,] > 0.05 & fa[i,] < 0.9
        if (!sum(idx) >= min.normals) return(c(0, 0, 0, 0))
        c(sum(ref[i,idx]), sum(alt[i,idx]), sum(idx), mean(fa[i, idx]))
    })
    # Add an average "normal" SNP (average coverage and allelic fraction > 0.4)
    # as empirical prior
    psMappingBias <- .adjustEmpBayes(psMappingBias) * 2
    ponCntHits <- apply(geno(nvcf)$AD, 1, function(x)
        sum(!is.na(unlist(x))) / 2)

    tmp <- rowRanges(nvcf)
    mcols(tmp) <- NULL
    tmp$bias <- psMappingBias
    tmp$pon.count <- ponCntHits
    tmp
}
    
.smoothVectorByChromosome <- function(x, chr, smooth.n) {
    .filter <- function(x, ...) {
        if (length(x) < smooth.n * 5) return(x)
        stats::filter(x, ...)
    }
    fN <- rep(1 / smooth.n, smooth.n)
    y <- do.call(c, lapply(
            split(x, factor(as.character(chr), levels = unique(chr))),
            .filter, fN, sides = 2))

    y[is.na(y)] <- x[is.na(y)]
    as.numeric(y)
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
    xg <- x[, x[4, ] > 0.4]
    # calculate the average number of ref and alt reads per sample
    shape1 <- sum(xg[1, ]) / sum(xg[3, ])
    shape2 <- sum(xg[2, ]) / sum(xg[3, ])
    # add those as empirical bayes estimate to all SNPs
    x[1, ] <- x[1, ] + shape1
    x[2, ] <- x[2, ] + shape2
    # get the alt allelic fraction for all SNPs
    apply(x, 2, function(y) y[2] / sum(y[1:2]))
}
