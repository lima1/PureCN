\name{segmentationCBS}
\alias{segmentationCBS}
\title{CBS segmentation}
\description{The default segmentation function. This function is called via the 
fun.segmentation argument of runAbsoluteCN. The arguments are passed
via args.segmentation.}
\usage{segmentationCBS(normal, tumor, log.ratio, plot.cnv, 
    coverage.cutoff, sampleid = sampleid, exon.weight.file = NULL, 
    alpha = 0.005, vcf = NULL, tumor.id.in.vcf = 1, 
    verbose = TRUE, ...)}
\arguments{
  \item{normal}{GATK coverage file for normal sample.}
  \item{tumor}{GATK coverage file for tumor sample.}
  \item{log.ratio}{Copy number log-ratios, one for each exon in coverage file.}
  \item{plot.cnv}{Segmentation plots.}
  \item{coverage.cutoff}{Minimum coverage in both normal and tumor.}
  \item{sampleid}{Sample id, used in output files.}
  \item{exon.weight.file}{Can be used to assign weights to exons.}
  \item{alpha}{Alpha value for CBS, see documentation for the segment function.}
  \item{vcf}{Optional VCF object with germline allelic ratios.}
  \item{tumor.id.in.vcf}{Id of tumor in case multiple samples are stored in VCF.}
  \item{verbose}{Verbose output.}
  \item{\dots}{Additional parameters passed to a modified version of CNV.analyze function
of the ExomeCNV package.}
}

\value{A list with elements
\item{seg}{The segmentation.}
\item{size}{The size of all segments (in base pairs).}}

\author{Markus Riester}




\examples{
gatk.normal.file <- system.file("extdata", "example_normal.txt", package="PureCN")
gatk.tumor.file <- system.file("extdata", "example_tumor.txt", package="PureCN")
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", package="PureCN")

# speed-up the runAbsoluteCN by using the stored grid-search.
# (purecn.example.output$candidates).
data(purecn.example.output)

# The max.candidate.solutions is set to a very low value only to speed-up this example.
# This is not a good idea for real samples.
ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, gatk.tumor.file=gatk.tumor.file, 
   vcf.file=vcf.file, sampleid='Sample1', gc.gene.file=gc.gene.file, 
   candidates=purecn.example.output$candidates, max.candidate.solutions=2,
   fun.segmentation=segmentationCBS, args.segmentation=list(alpha=0.001))
}
