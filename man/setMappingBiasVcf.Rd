\name{setMappingBiasVcf}
\alias{setMappingBiasVcf}
\title{Set Mapping Bias VCF}
\description{Function to set mapping  bias for each
variant in the provided \code{CollapsedVCF} object.
Currently, it returns the same value for all variants.}
\usage{setMappingBiasVcf(vcf, tumor.id.in.vcf = NULL, max.bias = 1, 
    normal.panel.vcf.file = NULL, min.normals = 10, 
    verbose = TRUE)}
\arguments{
  \item{vcf}{\code{CollapsedVCF} object, read in with the \code{readVcf} function 
from the VariantAnnotation package.}
  \item{tumor.id.in.vcf}{Id of tumor in case multiple samples are stored in VCF.}
  \item{max.bias}{Defines the maximum value for the mapping bias scaling factor.
The default of 1 assumes that the reference allele can never have
a lower mappability than the alt allele.}
  \item{normal.panel.vcf.file}{Combined VCF file of a panel of normals, expects allelic fractions
as FA genotype field. Should be compressed and indexed with bgzip and 
tabix, respectively.}
  \item{min.normals}{Minimum number of normals with heterozygous SNP. }
  \item{verbose}{Verbose output.}
}

\value{A \code{numeric(nrow(vcf))} vector with the mapping bias of 
for each variant in the \code{CollapsedVCF}. Mapping bias is expected as
scaling factor. Adjusted allelic fraction is 
(observed allelic fraction)/(mapping bias).   }

\author{Markus Riester}




\examples{
# This function is typically only called by runAbsoluteCN via the 
# fun.setMappingBiasVcf and args.setMappingBiasVcf comments.
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
vcf <- readVcf(vcf.file, "hg19")
vcf.bias <- setMappingBiasVcf(vcf)        
}
