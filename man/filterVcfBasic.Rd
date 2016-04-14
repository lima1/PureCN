\name{filterVcfBasic}
\alias{filterVcfBasic}
\title{filterVcfBasic}
\description{Function to remove artifacts and low confidence/quality calls from VCF file.}
\usage{filterVcfBasic(vcf, tumor.id.in.vcf = NULL, use.somatic.status = TRUE, 
    snp.blacklist = NULL, af.range = c(0.03, 0.97), 
    contamination.cutoff = c(0.05, 0.075), coverage.cutoff = 20, 
    min.supporting.reads = 3, verbose = TRUE)}
\arguments{
  \item{vcf}{VCF object, read in with the readVcf function from the VariantAnnotation package.}
  \item{tumor.id.in.vcf}{The tumor id in the VCF file, optional.}
  \item{use.somatic.status}{If somatic status and germline data is available, then clean up germline calls. }
  \item{snp.blacklist}{CSV file with SNP ids with expected allelic fraction significantly different from 0.5 in diploid genomes. Can be an array of lists. Can be created with createSNPBlacklist.}
  \item{af.range}{Exclude SNPs with allelic fraction smaller or greater than the two values, respectively.}
  \item{contamination.cutoff}{Count SNPs in dbSNP with allelic fraction smaller than the first value, if found on most chromosomes, remove all with AF smaller than the second value.}
  \item{coverage.cutoff}{Minimum coverage in tumor. Variants with lower coverage are ignored.}
  \item{min.supporting.reads}{Minimum number of reads supporting the alt allele.}
  \item{verbose}{
}
}

\value{A list with elements vcf, flag and flag_comment. "vcf" contains the filtered 
CollapsedVCF, "flag" a flag if problems were identified, further described in
"flag_comment".    }

\author{Markus Riester}




\examples{
### This function is typically only called by runAbsolute via the fun.filterVcf and 
### args.filterVcf comments.
library(VariantAnnotation)    
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
vcf <- readVcf(vcf.file, "hg19")
vcf.filtered <- filterVcfBasic(vcf)        
}
