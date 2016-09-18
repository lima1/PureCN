\name{setMappingBiasVcf}
\alias{setMappingBiasVcf}
\title{Set Mapping Bias VCF}
\description{Function to set mapping  bias for each
variant in the provided \code{CollapsedVCF} object.
Currently, it returns the same value for all variants.}
\usage{setMappingBiasVcf(vcf, tumor.id.in.vcf = NULL, verbose = TRUE)}
\arguments{
  \item{vcf}{\code{CollapsedVCF} object, read in with the \code{readVcf} function 
from the VariantAnnotation package.}
  \item{tumor.id.in.vcf}{Id of tumor in case multiple samples are stored in VCF.}
  \item{verbose}{Verbose output.}
}

\value{A \code{numeric(nrow(vcf))} vector with the mapping bias of 
for each variant in the \code{CollapsedVCF}.}

\author{Markus Riester}




\examples{
# This function is typically only called by runAbsoluteCN via the 
# fun.setPriorVcf and args.setPriorVcf comments.
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
vcf <- readVcf(vcf.file, "hg19")
vcf.bias <- setMappingBiasVcf(vcf)        
}
