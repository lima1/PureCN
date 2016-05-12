\name{getSexFromVcf}
\alias{getSexFromVcf}
\title{Get sample sex from VCF file}
\description{This function detects non-random distribution of homozygous
variants on chromosome X compared to all other chromosomes.}
\usage{getSexFromVcf(vcf, tumor.id.in.vcf = NULL, min.or = 4, 
    min.or.na = 2, max.pv = 1e-04, verbose = TRUE)}
\arguments{
  \item{vcf}{CollapsedVCF object, read in with the readVcf function 
from the VariantAnnotation package.}
  \item{tumor.id.in.vcf}{The tumor id in the CollapsedVCF (optional).}
  \item{min.or}{Minimum odds-ratio to call sample as male.}
  \item{min.or.na}{Minimum odds-ratio to not call a sample.}
  \item{max.pv}{Maximum Fisher's exact p-value to call sample as male.}
  \item{verbose}{Verbose output.}
}

\value{Returns "M" for male, "F" for female, or NA if unknown.    }

\author{Markus Riester}




\examples{
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
vcf <- readVcf(vcf.file, "hg19")
# This example vcf is already filtered and contains no homozygous calls,
# which are necessary for determining sex from chromosome X.
getSexFromVcf(vcf)
}
