\name{filterVcfBasic}
\alias{filterVcfBasic}
\title{Basic VCF filter function}
\description{Function to remove artifacts and low confidence/quality 
variant calls. }
\usage{filterVcfBasic(vcf, tumor.id.in.vcf = NULL, use.somatic.status = TRUE, 
    snp.blacklist = NULL, af.range = c(0.03, 0.97), 
    contamination.cutoff = c(0.05, 0.075), min.coverage = 15, 
    min.base.quality = 25, min.supporting.reads = NULL, 
    error = 0.001, target.granges = NULL, remove.off.target.snvs = TRUE, 
    interval.padding = 50, verbose = TRUE)}
\arguments{
  \item{vcf}{\code{CollapsedVCF} object, read in with the \code{readVcf} 
function from the VariantAnnotation package.}
  \item{tumor.id.in.vcf}{The tumor id in the \code{CollapsedVCF} (optional).}
  \item{use.somatic.status}{If somatic status and germline data is available, then use this 
information to remove non-heterozygous germline SNPs or germline SNPS 
with biased allelic fractions.}
  \item{snp.blacklist}{CSV file with SNP ids with expected allelic fraction 
significantly different from 0.5 in diploid genomes. Can be an array of 
lists. The function \code{\link{createSNPBlacklist}} can provide 
appropriate black lists. Can also be a BED file (either tab or 
comma separated) of  blacklisted genomic regions (columns 1-3: 
chromosome, start, end).}
  \item{af.range}{Exclude SNPs with allelic fraction smaller or greater than the 
two values, respectively. The higher value removes homozygous SNPs,
which potentially have allelic fractions smaller than 1 due to artifacts
or contamination. If a matched normal is available, this value is ignored,
because homozygosity can be confirmed in the normal.}
  \item{contamination.cutoff}{Count SNPs in dbSNP with allelic fraction smaller than the 
first value, if found on most chromosomes, remove all with AF smaller than
the second value.}
  \item{min.coverage}{Minimum coverage in tumor. Variants with lower coverage are ignored.}
  \item{min.base.quality}{Minimim base quality in tumor. Requires a \code{BQ} genotype field 
in the VCF.}
  \item{min.supporting.reads}{Minimum number of reads supporting the alt allele. 
If \code{NULL}, calculate based on coverage and assuming sequencing error 
of 10^-3.}
  \item{error}{Estimated sequencing error rate. Used to calculate minimum
number of supporting reads using \code{\link{calculatePowerDetectSomatic}}.}
  \item{target.granges}{\code{GenomicRanges} object specifiying the target postions.
Used to remove off-target reads. If \code{NULL}, do not check whether
variants are on or off-target.}
  \item{remove.off.target.snvs}{If set to a true value, will remove all SNVs outside the 
covered regions.}
  \item{interval.padding}{Include variants in the interval flanking regions of the 
specified size in bp. Requires \code{target.granges}.}
  \item{verbose}{Verbose output.}
}

\value{A list with elements
\item{vcf}{The filtered \code{CollapsedVCF} object.}
\item{flag}{A flag (\code{logical(1)}) if problems were identified.}
\item{flag_comment}{A comment describing the flagging.}}

\author{Markus Riester}



\seealso{\code{\link{calculatePowerDetectSomatic}}}
\examples{
# This function is typically only called by runAbsolute via the 
# fun.filterVcf and args.filterVcf comments.
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
vcf <- readVcf(vcf.file, "hg19")
vcf.filtered <- filterVcfBasic(vcf)        
}
