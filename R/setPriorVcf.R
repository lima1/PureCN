#' Set Somatic Prior VCF
#' 
#' Function to set prior for somatic mutation status for each variant in the
#' provided \code{CollapsedVCF} object.
#' 
#' 
#' @param vcf \code{CollapsedVCF} object, read in with the \code{readVcf}
#' function from the VariantAnnotation package.
#' @param prior.somatic Prior probabilities for somatic mutations. First value
#' is for the case when no matched normals are available and the variant is not
#' in dbSNP (second value). Third value is for variants with MuTect somatic
#' call. Different from 1, because somatic mutations in segments of copy number
#' 0 have 0 probability and artifacts can thus have dramatic influence on
#' likelihood score. Forth value is for variants not labeled as somatic by
#' MuTect. Last two values are optional, if vcf contains a flag Cosmic.CNT, it
#' will set the prior probability for variants with CNT > 2 to the first of
#' those values in case of no matched normal available (0.995 default).  Final
#' value is for the case that variant is in both dbSNP and COSMIC > 2.
#' @param tumor.id.in.vcf Id of tumor in case multiple samples are stored in
#' VCF.
#' @param min.cosmic.cnt Minimum number of hits in the COSMIC database to 
#' call variant as likely somatic.
#' @param DB.info.flag Flag in INFO of VCF that marks presence in common
#' germline databases. Defaults to \code{DB} that may contain somatic variants
#' if it is from an unfiltered dbSNP VCF.
#' @return A \code{numeric(nrow(vcf))} vector with the prior probability of
#' somatic status for each variant in the \code{CollapsedVCF}.
#' @author Markus Riester
#' @examples
#' 
#' # This function is typically only called by runAbsoluteCN via the 
#' # fun.setPriorVcf and args.setPriorVcf comments.
#' vcf.file <- system.file("extdata", "example_vcf.vcf.gz", package="PureCN")
#' vcf <- readVcf(vcf.file, "hg19")
#' vcf.priorsomatic <- setPriorVcf(vcf)        
#' 
#' @export setPriorVcf
setPriorVcf <- function(vcf, prior.somatic = c(0.5, 0.0005, 0.999, 0.0001,
                                               0.995, 0.5),
                        tumor.id.in.vcf = NULL, min.cosmic.cnt = 4, 
                        DB.info.flag = "DB") {
    if (is.null(tumor.id.in.vcf)) {
        tumor.id.in.vcf <- .getTumorIdInVcf(vcf)
    }
    if (!is.null(info(vcf)$SOMATIC)) {
         tmp <- prior.somatic
         prior.somatic <- ifelse(info(vcf)$SOMATIC,
            prior.somatic[3],prior.somatic[4])

         flog.info("Found SOMATIC annotation in VCF.")
         flog.info("Setting somatic prior probabilities for somatic variants to %f or to %f otherwise.", 
            tmp[3], tmp[4])
    } else {
         tmp <- prior.somatic
         prior.somatic <- ifelse(info(vcf)[[DB.info.flag]],
            prior.somatic[2], prior.somatic[1])
         if (!is.null(info(vcf)$Cosmic.CNT)) {
             flog.info("Found COSMIC annotation in VCF.")
             flog.info("Setting somatic prior probabilities for hits to %f or to %f if in both COSMIC and dbSNP.", 
                tmp[5], tmp[6])

             prior.somatic[which(info(vcf)$Cosmic.CNT >= min.cosmic.cnt)] <- tmp[5]
             prior.somatic[which(info(vcf)$Cosmic.CNT >= min.cosmic.cnt & 
                info(vcf)[[DB.info.flag]])] <- tmp[6]
         } else {
             flog.info("Setting somatic prior probabilities for dbSNP hits to %f or to %f otherwise.", 
                tmp[2], tmp[1])
         }      
    }     
    prior.somatic
}
