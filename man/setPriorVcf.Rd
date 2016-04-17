\name{setPriorVcf}
\alias{setPriorVcf}
\title{Set Somatic Prior VCF}
\description{Function to set prior for somatic mutation status for each
variant in the provided CollapsedVCF object.}
\usage{setPriorVcf(vcf, prior.somatic = c(0.5, 5e-04, 0.999, 
    1e-04, 0.95, 0.01), verbose = TRUE)}
\arguments{
  \item{vcf}{VCF object, read in with the readVcf function from the 
VariantAnnotation package.}
  \item{prior.somatic}{Prior probabilities for somatic mutations. First value is for 
the case when no matched normals are available and the variant is not in 
dbSNP (second value). Third value is for variants with Mutect somatic call.
Different from 1, because somatic mutations in segments of copy number 0 
have 0 probability and artifacts can thus have dramatic influence on 
likelihood score. Forth value is for variants not labeled as somatic by 
Mutect. Last two values are optional, if vcf contains a flag Cosmic.CNT, 
it will set the prior probability for variants with CNT > 2 to the first 
of those values in case of no matched normal available (0.95 default). 
Final value is for the case that variant is in both dbSNP and COSMIC > 2. }
  \item{verbose}{Verbose output.}
}

\value{A vector with the prior probability of somatic status for each 
variant in the CollapsedVCF.}

\author{Markus Riester}




\examples{
# This function is typically only called by runAbsoluteCN via the 
# fun.setPriorVcf and args.setPriorVcf comments.
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
vcf <- readVcf(vcf.file, "hg19")
vcf.priorsomatic <- setPriorVcf(vcf)        
}
