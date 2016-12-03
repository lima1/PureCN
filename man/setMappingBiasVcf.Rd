\name{setMappingBiasVcf}
\alias{setMappingBiasVcf}
\title{Set Mapping Bias VCF}
\description{Function to set mapping  bias for each
variant in the provided \code{CollapsedVCF} object.
By default, it returns the same value for all variants, but a
pool of normal samples can be provided for position-specific
mapping bias calculation.}
\usage{setMappingBiasVcf(vcf, tumor.id.in.vcf = NULL, normal.panel.vcf.file = NULL, 
    min.normals = 5, smooth = TRUE, smooth.n = 5, verbose = TRUE)}
\arguments{
  \item{vcf}{\code{CollapsedVCF} object, read in with the \code{readVcf} function 
from the VariantAnnotation package.}
  \item{tumor.id.in.vcf}{Id of tumor in case multiple samples are stored in VCF.}
  \item{normal.panel.vcf.file}{Combined VCF file of a panel of normals, expects allelic fractions
as FA genotype field. Should be compressed and indexed with bgzip and 
tabix, respectively.}
  \item{min.normals}{Minimum number of normals with heterozygous SNP for calculating
position-specific mapping bias. Requires \code{normal.panel.vcf.file}.  }
  \item{smooth}{Impute mapping bias of variants not found in the panel by smoothing
of neighboring SNPs. Requires \code{normal.panel.vcf.file}.  }
  \item{smooth.n}{Number of neighboring variants used for smoothing.}
  \item{verbose}{Verbose output.}
}

\value{A \code{numeric(nrow(vcf))} vector with the mapping bias of 
for each variant in the \code{CollapsedVCF}. Mapping bias is expected as
scaling factor. Adjusted allelic fraction is 
(observed allelic fraction)/(mapping bias). Maximum scaling factor is 1
and means no bias. }

\author{Markus Riester}




\examples{
# This function is typically only called by runAbsoluteCN via 
# fun.setMappingBiasVcf and args.setMappingBiasVcf.
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
vcf <- readVcf(vcf.file, "hg19")
vcf.bias <- setMappingBiasVcf(vcf)        
}
