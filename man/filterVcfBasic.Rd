\name{filterVcfBasic}
\alias{filterVcfBasic}
\title{Basic VCF filter function}
\description{Function to remove artifacts and low confidence/quality 
variant calls.}
\usage{filterVcfBasic(vcf, tumor.id.in.vcf = NULL, use.somatic.status = TRUE, 
    snp.blacklist = NULL, af.range = c(0.03, 0.97), 
    contamination.cutoff = c(0.05, 0.075), coverage.cutoff = 20, 
    min.supporting.reads = NULL, error = 0.001, verbose = TRUE)}
\arguments{
  \item{vcf}{CollapsedVCF object, read in with the readVcf function 
from the VariantAnnotation package.}
  \item{tumor.id.in.vcf}{The tumor id in the CollapsedVCF (optional).}
  \item{use.somatic.status}{If somatic status and germline data is available, then use this 
information to remove non-heterozygous germline SNPs or germline SNPS 
with biased allelic fractions.}
  \item{snp.blacklist}{CSV file with SNP ids with expected allelic fraction 
significantly different from 0.5 in diploid genomes. Can be an array of 
lists. The function createSNPBlacklist can provide appropriate black 
lists. Can also be a BED file (either tab or comma separated) of 
blacklisted genomic regions (columns 1-3: chromosome, start, end).}
  \item{af.range}{Exclude SNPs with allelic fraction smaller or greater than the 
two values, respectively.}
  \item{contamination.cutoff}{Count SNPs in dbSNP with allelic fraction smaller than the 
first value, if found on most chromosomes, remove all with AF smaller than
the second value.}
  \item{coverage.cutoff}{Minimum coverage in tumor. Variants with lower coverage are ignored.}
  \item{min.supporting.reads}{Minimum number of reads supporting the alt allele. 
If NULL, calculate based on coverage and assuming sequencing error 
of 10^-3.}
  \item{error}{Estimated sequencing error rate. Used to calculate minimum
number of supporting reads using calculatePowerDetectSomatic.}
  \item{verbose}{
}
}

\value{A list with elements
\item{vcf}{The filtered CollapsedVCF object.}
\item{flag}{A flag (TRUE/FALSE) if problems were identified.}
\item{flag_comment}{A comment describing the flagging.}}

\author{Markus Riester}




\examples{
# This function is typically only called by runAbsolute via the 
# fun.filterVcf and args.filterVcf comments.
library(VariantAnnotation)    
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
vcf <- readVcf(vcf.file, "hg19")
vcf.filtered <- filterVcfBasic(vcf)        
}
