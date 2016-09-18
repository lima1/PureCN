setMappingBiasVcf <- structure(function(# Set Mapping Bias VCF
### Function to set mapping  bias for each
### variant in the provided \code{CollapsedVCF} object.
### Currently, it returns the same value for all variants.
vcf,
### \code{CollapsedVCF} object, read in with the \code{readVcf} function 
### from the VariantAnnotation package.
tumor.id.in.vcf=NULL,
### Id of tumor in case multiple samples are stored in VCF.
verbose=TRUE
### Verbose output.
) {
    if (is.null(tumor.id.in.vcf)) {
        tumor.id.in.vcf <- .getTumorIdInVcf(vcf)
    }
    mappingBias <- 1
    if (!is.null(info(vcf)$SOMATIC) && ncol(vcf)>1) {
         normal.id.in.vcf <- .getNormalIdInVcf(vcf, tumor.id.in.vcf)
         faAll <- as.numeric(geno(vcf)$FA[,normal.id.in.vcf])
         mappingBias <- mean(faAll, na.rm=TRUE)*2
         if (verbose) message("Found SOMATIC annotation in VCF. ",
            "Setting mapping bias to ", round(mappingBias, digits=3)) 
    }     
    rep(mappingBias, nrow(vcf))
### A \code{numeric(nrow(vcf))} vector with the mapping bias of 
### for each variant in the \code{CollapsedVCF}.
},ex=function() {
# This function is typically only called by runAbsoluteCN via the 
# fun.setPriorVcf and args.setPriorVcf comments.
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
vcf <- readVcf(vcf.file, "hg19")
vcf.bias <- setMappingBiasVcf(vcf)        
})  
