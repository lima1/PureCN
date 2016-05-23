\name{getSexFromVcf}
\alias{getSexFromVcf}
\title{Get sample sex from a VCF file}
\description{This function detects non-random distribution of homozygous
variants on chromosome X compared to all other chromosomes.
A non-significant Fisher's exact p-value indicates more than one
chromosome X copy. This function is called in runAbsoluteCN as 
sanity check when a VCF is provided. It is also useful for determining
sex when no sex marker genes on chrY (e.g. AMELY) are available.}
\usage{getSexFromVcf(vcf, tumor.id.in.vcf = NULL, min.or = 4, 
    min.or.na = 2.5, max.pv = 0.001, homozygous.cutoff = 0.95, 
    af.cutoff = 0.03, verbose = TRUE)}
\arguments{
  \item{vcf}{CollapsedVCF object, read in with the readVcf function 
from the VariantAnnotation package.}
  \item{tumor.id.in.vcf}{The tumor id in the CollapsedVCF (optional).}
  \item{min.or}{Minimum odds-ratio to call sample as male. If p-value is
not significant due to a small number of SNPs on chromosome X,
sample will be called as NA even when odds-ratio exceeds this cutoff.}
  \item{min.or.na}{Minimum odds-ratio to not call a sample. Odds-ratios in the
range min.or.na to min.or define a grey area in which samples
are not called. Contamination can be a source of ambiguous calls.}
  \item{max.pv}{Maximum Fisher's exact p-value to call sample as male.}
  \item{homozygous.cutoff}{Minimum allelic fraction to call position homozygous.}
  \item{af.cutoff}{Remove all SNVs with allelic fraction lower than the
specified value. }
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
