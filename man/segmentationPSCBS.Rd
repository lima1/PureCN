\name{segmentationPSCBS}
\alias{segmentationPSCBS}
\title{PSCBS segmentation}
\description{Alternative segmentation function using the \code{PSCBS} package. 
This function is called via the 
\code{fun.segmentation} argument of \code{\link{runAbsoluteCN}}. 
The arguments are passed via \code{args.segmentation}.}
\usage{segmentationPSCBS(normal, tumor, log.ratio, seg, plot.cnv, 
    min.coverage, sampleid = sampleid, exon.weight.file = NULL, 
    flavor = "tcn&dh", tauA = 0.03, vcf = NULL, tumor.id.in.vcf = 1, 
    normal.id.in.vcf = NULL, max.segments = NULL, prune.hclust.h = NULL, 
    prune.hclust.method = "ward.D", chr.hash = NULL, 
    verbose = TRUE, ...)}
\arguments{
  \item{normal}{GATK coverage data for normal sample.}
  \item{tumor}{GATK coverage data for tumor sample.}
  \item{log.ratio}{Copy number log-ratios, one for each exon in coverage file.}
  \item{seg}{If segmentation was provided by the user, this data structure will contain
this segmentation. Useful for minimal segmentation functions. Otherwise
PureCN will re-segment the data. This segmentation function ignores this
user provided segmentation.}
  \item{plot.cnv}{Segmentation plots.}
  \item{min.coverage}{Minimum coverage in both normal and tumor.}
  \item{sampleid}{Sample id, used in output files.}
  \item{exon.weight.file}{Can be used to assign weights to exons. NOT SUPPORTED YET.}
  \item{flavor}{Flavor value for PSCBS. See segmentByNonPairedPSCBS.}
  \item{tauA}{tauA argument for PSCBS. See segmentByNonPairedPSCBS.}
  \item{vcf}{Optional VCF object with germline allelic ratios.}
  \item{tumor.id.in.vcf}{Id of tumor in case multiple samples are stored in VCF.}
  \item{normal.id.in.vcf}{Id of normal in in VCF. If NULL, use unpaired PSCBS.}
  \item{max.segments}{If not NULL, try a higher undo.SD parameter if number of
segments exceeds the threshold.}
  \item{prune.hclust.h}{Height in the \code{hclust} pruning step. Increasing this value 
will merge segments more aggressively. If NULL, try to find a sensible
default.}
  \item{prune.hclust.method}{Cluster method used in the \code{hclust} pruning step. See
documentation for the \code{hclust} function.}
  \item{chr.hash}{Mapping of non-numerical chromsome names to numerical names
(e.g. chr1 to 1, chr2 to 2, etc.). If NULL, assume chromsomes
are properly ordered.   }
  \item{verbose}{Verbose output.}
  \item{\dots}{Additional parameters passed to the segmentByNonPairedPSCBS function.}
}

\value{A list with elements seg and size. "seg" contains the 
segmentation, "size" the size of all segments in base pairs.    }

\author{Markus Riester}



\seealso{\code{\link{runAbsoluteCN}}}
\examples{
normal.coverage.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")
vcf.file <- system.file("extdata", "example_vcf.vcf", 
    package="PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
    package="PureCN")

# The max.candidate.solutions, max.ploidy and test.purity parameters are set to
# non-default values to speed-up this example.  This is not a good idea for real
# samples.
 ret <-runAbsoluteCN(normal.coverage.file=normal.coverage.file, 
     tumor.coverage.file=tumor.coverage.file, vcf.file=vcf.file, genome="hg19",
     sampleid='Sample1', gc.gene.file=gc.gene.file,
     fun.segmentation=segmentationPSCBS, max.ploidy=4,
     test.purity=seq(0.3,0.7,by=0.05), max.candidate.solutions=1)
}
