setMappingBiasVcf <- structure(function(# Set Mapping Bias VCF
### Function to set mapping  bias for each
### variant in the provided \code{CollapsedVCF} object.
### By default, it returns the same value for all variants, but a
### pool of normal samples can be provided for position-specific
### mapping bias calculation.
vcf,
### \code{CollapsedVCF} object, read in with the \code{readVcf} function 
### from the VariantAnnotation package.
tumor.id.in.vcf=NULL,
### Id of tumor in case multiple samples are stored in VCF.
normal.panel.vcf.file=NULL,
### Combined VCF file of a panel of normals, expects allelic fractions
### as FA genotype field. Should be compressed and indexed with bgzip and 
### tabix, respectively.
min.normals=5,
### Minimum number of normals with heterozygous SNP for calculating
### position-specific mapping bias. Requires \code{normal.panel.vcf.file}.  
smooth=TRUE,
### Impute mapping bias of variants not found in the panel by smoothing
### of neighboring SNPs. Requires \code{normal.panel.vcf.file}.  
smooth.n=5,
### Number of neighboring variants used for smoothing.
verbose=TRUE
### Verbose output.
) {
    if (is.null(tumor.id.in.vcf)) {
        tumor.id.in.vcf <- .getTumorIdInVcf(vcf)
    }
    mappingBias <- 1
    if (!is.null(info(vcf)$SOMATIC) && ncol(vcf)>1) {
         normal.id.in.vcf <- .getNormalIdInVcf(vcf, tumor.id.in.vcf)
         faAll <- as.numeric(geno(vcf)$FA[!info(vcf)$SOMATIC,normal.id.in.vcf])
         mappingBias <- mean(faAll, na.rm=TRUE)*2
         if (verbose) message("Found SOMATIC annotation in VCF. ",
            "Setting mapping bias to ", round(mappingBias, digits=3)) 
    }     
    if (is.null(info(vcf)$SOMATIC) && is.null(normal.panel.vcf.file)) {
        message("VCF does not contain somatic status and no ",
        "normal.panel.vcf.file provided.")
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
    nvcf <- .readNormalPanelVcfLarge(vcf, normal.panel.vcf.file, 
        verbose=verbose)
    if (nrow(nvcf) < 1) {
        warning("setMappingBiasVcf: no hits in ", normal.panel.vcf.file, ".")
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
### A \code{numeric(nrow(vcf))} vector with the mapping bias of 
### for each variant in the \code{CollapsedVCF}. Mapping bias is expected as
### scaling factor. Adjusted allelic fraction is 
### (observed allelic fraction)/(mapping bias). Maximum scaling factor is 1
### and means no bias. 
},ex=function() {
# This function is typically only called by runAbsoluteCN via 
# fun.setMappingBiasVcf and args.setMappingBiasVcf.
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
vcf <- readVcf(vcf.file, "hg19")
vcf.bias <- setMappingBiasVcf(vcf)        
}) 

.smoothVectorByChromosome <- function(x, chr, smooth.n) {
    .filter <- function(x, ...) {
        if (length(x) < smooth.n*5) return(x)
        filter(x, ...)
    }       
    fN <- rep(1/smooth.n, smooth.n)
    y <- do.call(c, lapply(split(x, factor(as.character(chr), levels=unique(chr))), .filter, fN, sides=2))
    y[is.na(y)] <- x[is.na(y)]
    as.numeric(y)
}

.readNormalPanelVcfLarge <- function(vcf, normal.panel.vcf.file, max.file.size=1, verbose) {
    genome <- genome(vcf)[1]    
    if (file.size(normal.panel.vcf.file)/1000^3 > max.file.size || nrow(vcf)< 1000) {
        if (verbose) message("Scanning ", normal.panel.vcf.file, "...")
        nvcf <- readVcf(TabixFile(normal.panel.vcf.file), genome=genome, 
            ScanVcfParam(which = rowRanges(vcf), info=NA, fixed=NA, 
            geno="FA"))
    } else {
        if (verbose) message("Loading ", normal.panel.vcf.file, "...")
        nvcf <- readVcf(normal.panel.vcf.file, genome=genome,
            ScanVcfParam(info=NA, fixed=NA, geno="FA"))
        nvcf <- subsetByOverlaps(nvcf, rowRanges(vcf))
    }    
    nvcf
}    
