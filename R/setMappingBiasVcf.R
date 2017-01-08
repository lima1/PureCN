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
#' expects allelic fractions as FA genotype field. Should be compressed and
#' indexed with bgzip and tabix, respectively.
#' @param min.normals Minimum number of normals with heterozygous SNP for
#' calculating position-specific mapping bias. Requires
#' \code{normal.panel.vcf.file}.
#' @param smooth Impute mapping bias of variants not found in the panel by
#' smoothing of neighboring SNPs. Requires \code{normal.panel.vcf.file}.
#' @param smooth.n Number of neighboring variants used for smoothing.
#' @return A \code{numeric(nrow(vcf))} vector with the mapping bias of for each
#' variant in the \code{CollapsedVCF}. Mapping bias is expected as scaling
#' factor. Adjusted allelic fraction is (observed allelic fraction)/(mapping
#' bias). Maximum scaling factor is 1 and means no bias.
#' @author Markus Riester
#' @examples
#' 
#' # This function is typically only called by runAbsoluteCN via 
#' # fun.setMappingBiasVcf and args.setMappingBiasVcf.
#' vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
#' vcf <- readVcf(vcf.file, "hg19")
#' vcf.bias <- setMappingBiasVcf(vcf)        
#' 
#' @export setMappingBiasVcf
setMappingBiasVcf <- function(vcf, tumor.id.in.vcf = NULL,
normal.panel.vcf.file = NULL, min.normals = 5, smooth = TRUE, smooth.n = 5) {

    if (is.null(tumor.id.in.vcf)) {
        tumor.id.in.vcf <- .getTumorIdInVcf(vcf)
    }
    mappingBias <- 1
    if (!is.null(info(vcf)$SOMATIC) && ncol(vcf)>1) {
         normal.id.in.vcf <- .getNormalIdInVcf(vcf, tumor.id.in.vcf)
         faAll <- as.numeric(geno(vcf)$FA[!info(vcf)$SOMATIC,normal.id.in.vcf])
         mappingBias <- mean(faAll, na.rm=TRUE)*2
         flog.info("Found SOMATIC annotation in VCF. Setting mapping bias to %.3f.", 
            mappingBias) 
    }     
    if (is.null(info(vcf)$SOMATIC) && is.null(normal.panel.vcf.file)) {
        flog.info(
            "VCF does not contain somatic status. For best results, consider%s%s",
            "providing normal.panel.vcf.file when matched normals are not ",
            "available.")
    }    
    tmp <- rep(mappingBias, nrow(vcf)) 
    # Defines the maximum value for the mapping bias scaling factor.
    # 1 assumes that the reference allele can never have
    # a lower mappability than the alt allele.
    max.bias <- 1
    tmp[tmp>max.bias] <- max.bias
    if (is.null(normal.panel.vcf.file)) {
        return(tmp)
    } 
    nvcf <- .readNormalPanelVcfLarge(vcf, normal.panel.vcf.file)
    if (nrow(nvcf) < 1) {
        flog.warn("setMappingBiasVcf: no hits in %s.", normal.panel.vcf.file)
        return(tmp)
    }

    psMappingBias <- apply(geno(nvcf)$FA, 1, function(x) { 
        x <- unlist(x)
        x <- subset(x, !is.na(x) & x>0.1 &x < 0.9)
        ifelse(length(x) >= min.normals, mean(x),0.5)
    })*2
    
    ov <- findOverlaps(vcf, nvcf, select="first")

    tmp[!is.na(ov)] <- psMappingBias[ov][!is.na(ov)]
    tmp[tmp>max.bias] <- max.bias
    if (smooth) {
        tmpSmoothed <- .smoothVectorByChromosome(tmp, as.character(seqnames(vcf)), 
            smooth.n)
        # only smooth variants not found in database
        tmp[is.na(ov)] <- tmpSmoothed[is.na(ov)]
    }
    return(tmp)
}

.smoothVectorByChromosome <- function(x, chr, smooth.n) {
    .filter <- function(x, ...) {
        if (length(x) < smooth.n*5) return(x)
        stats::filter(x, ...)
    }       
    fN <- rep(1/smooth.n, smooth.n)
    y <- do.call(c, lapply(split(x, factor(as.character(chr), levels=unique(chr))), .filter, fN, sides=2))
    y[is.na(y)] <- x[is.na(y)]
    as.numeric(y)
}

.readNormalPanelVcfLarge <- function(vcf, normal.panel.vcf.file, max.file.size=1) {
    genome <- genome(vcf)[1]    
    if (file.size(normal.panel.vcf.file)/1000^3 > max.file.size || nrow(vcf)< 1000) {
        flog.info("Scanning %s...", normal.panel.vcf.file)
        nvcf <- readVcf(TabixFile(normal.panel.vcf.file), genome=genome, 
            ScanVcfParam(which = rowRanges(vcf), info=NA, fixed=NA, 
            geno="FA"))
    } else {
        flog.info("Loading %s...", normal.panel.vcf.file)
        nvcf <- readVcf(normal.panel.vcf.file, genome=genome,
            ScanVcfParam(info=NA, fixed=NA, geno="FA"))
        nvcf <- subsetByOverlaps(nvcf, rowRanges(vcf))
    }    
    nvcf
}    
