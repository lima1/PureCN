#' Set Mapping Bias VCF
#' 
#' Function to set mapping bias for each variant in the provided
#' \code{CollapsedVCF} object. By default, it returns the same value for all
#' variants, but a mapping bias file can be provided for position-specific
#' mapping bias calculation.
#' 
#' 
#' @param vcf \code{CollapsedVCF} object, read in with the \code{readVcf}
#' function from the VariantAnnotation package.
#' @param tumor.id.in.vcf Id of tumor in case multiple samples are stored in
#' VCF.
#' @param mapping.bias.file A precomputed mapping bias database 
#' obtained by \code{\link{calculateMappingBiasVcf}}.
#' instead.
#' reference and alt counts as AD genotype field. Should be compressed and
#' @param smooth Impute mapping bias of variants not found in the panel by
#' smoothing of neighboring SNPs. Requires \code{mapping.bias.file}.
#' @param smooth.n Number of neighboring variants used for smoothing.
#' @return A \code{data.frame} with elements \item{bias}{A \code{numeric(nrow(vcf))} 
#' vector with the mapping bias of for each
#' variant in the \code{CollapsedVCF}. Mapping bias is expected as scaling
#' factor. Adjusted allelic fraction is (observed allelic fraction)/(mapping
#' bias). Maximum scaling factor is 1 and means no bias.}
#' \item{pon.count}{A \code{numeric(nrow(vcf))} vector with the number
#' of hits in the \code{mapping.bias.file}.}
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
#' @importFrom GenomeInfoDb genome<-
setMappingBiasVcf <- function(vcf, tumor.id.in.vcf = NULL, mapping.bias.file = NULL,
smooth = TRUE, smooth.n = 5) {
    if (is.null(tumor.id.in.vcf)) {
        tumor.id.in.vcf <- .getTumorIdInVcf(vcf)
    }

    mappingBias <- 1
    if (!is.null(info(vcf)$SOMATIC) && ncol(vcf) > 1) {
         normal.id.in.vcf <- .getNormalIdInVcf(vcf, tumor.id.in.vcf)
         faAll <- as.numeric(geno(vcf)$FA[!info(vcf)$SOMATIC, normal.id.in.vcf])
         mappingBias <- mean(faAll, na.rm=TRUE) * 2
         flog.info("Found SOMATIC annotation in VCF. Setting mapping bias to %.3f.",
            mappingBias)
    }
    if (is.null(info(vcf)$SOMATIC) && is.null(mapping.bias.file)) {
        flog.info(
            "VCF does not contain somatic status. For best results, consider%s%s",
            " providing mapping.bias.file when matched normals are not ",
            "available.")
    }
    tmp <- rep(mappingBias, nrow(vcf))
    # Defines the maximum value for the mapping bias scaling factor.
    # 1 assumes that the reference allele can never have
    # a lower mappability than the alt allele.
    max.bias <- 1.2
    tmp[tmp > max.bias] <- max.bias
    if (is.null(mapping.bias.file)) {
        return(data.frame(bias = tmp, shape1 = NA, shape2 = NA))
    }

    if (tolower(file_ext(mapping.bias.file)) == "rds") {
        flog.info("Loading mapping bias file %s...",
            basename(mapping.bias.file))
        mappingBias <- readRDS(mapping.bias.file)
    } else {
        .stopUserError("mapping.bias.file must be a file with *.rds suffix.")
    }
    .annotateMappingBias(tmp, vcf, mappingBias, max.bias, smooth, smooth.n)
}



.annotateMappingBias <- function(tmp, vcf, mappingBias, max.bias, smooth, smooth.n) {
    mappingBias <- .checkSeqlevelStyle(vcf, mappingBias, "mapping.bias.file", "vcf")
    .compareGenomes <- function(x, y) {
        gx <- genome(x)
        gy <- genome(y)
        isc <- intersect(names(gx), names(gy))
        identical(gx[isc], gy[isc])
     }       
    if (!.compareGenomes(vcf, mappingBias)) {
        genome(mappingBias) <- genome(vcf)
    }
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
    mu <- rep(NA, length(tmp))
    rho <- rep(NA, length(tmp))
    ov <- findOverlaps(vcf, mappingBias, select = "first")
    idx <- !is.na(ov)
    tmp[idx] <- mappingBias$bias[ov][idx]
    ponCnt[idx] <- mappingBias$pon.count[ov][idx]
    if (!is.null(mappingBias$mu)) {
        mu[idx] <- mappingBias$mu[ov][idx]
        rho[idx] <- mappingBias$rho[ov][idx]
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
                mu = mu, rho = rho))
}
    
