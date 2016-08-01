\name{segmentationCBS}
\alias{segmentationCBS}
\title{CBS segmentation}
\description{The default segmentation function. This function is called via the 
\code{fun.segmentation} argument of \code{\link{runAbsoluteCN}}. 
The arguments are passed via \code{args.segmentation}.}
\usage{segmentationCBS(normal, tumor, log.ratio, plot.cnv, 
    coverage.cutoff, sampleid = sampleid, exon.weight.file = NULL, 
    alpha = 0.005, undo.SD = NULL, vcf = NULL, tumor.id.in.vcf = 1, 
    normal.id.in.vcf = NULL, max.segments = NULL, prune.hclust.h = NULL, 
    prune.hclust.method = "ward.D", chr.hash = NULL, 
    verbose = TRUE)}
\arguments{
  \item{normal}{GATK coverage data for normal sample.}
  \item{tumor}{GATK coverage data for tumor sample.}
  \item{log.ratio}{Copy number log-ratios, one for each exon in coverage file.}
  \item{plot.cnv}{Segmentation plots.}
  \item{coverage.cutoff}{Minimum coverage in both normal and tumor.}
  \item{sampleid}{Sample id, used in output files.}
  \item{exon.weight.file}{Can be used to assign weights to exons.}
  \item{alpha}{Alpha value for CBS, see documentation for the \code{segment} function.}
  \item{undo.SD}{\code{undo.SD} for CBS, see documentation of the \code{segment} function.
If NULL, try to find a sensible default.}
  \item{vcf}{Optional \code{CollapsedVCF} object with germline allelic ratios.}
  \item{tumor.id.in.vcf}{Id of tumor in case multiple samples are stored in VCF.}
  \item{normal.id.in.vcf}{Id of normal in in VCF. Currently not used.}
  \item{max.segments}{If not \code{NULL}, try a higher \code{undo.SD} parameter 
if number of segments exceeds the threshold.}
  \item{prune.hclust.h}{Height in the \code{hclust} pruning step. Increasing this value 
will merge segments more aggressively. If NULL, try to find a sensible
default.}
  \item{prune.hclust.method}{Cluster method used in the \code{hclust} pruning step. See
documentation for the \code{hclust} function.}
  \item{chr.hash}{Mapping of non-numerical chromsome names to numerical names
(e.g. chr1 to 1, chr2 to 2, etc.). If \code{NULL}, assume chromsomes
are properly ordered.}
  \item{verbose}{Verbose output.}
}

\value{\code{data.frame} containing the segmentation.    }

\author{Markus Riester}



\seealso{\code{\link{runAbsoluteCN}}}
\examples{
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")
vcf.file <- system.file("extdata", "example_vcf.vcf", 
    package="PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
    package="PureCN")

# speed-up the runAbsoluteCN by using the stored grid-search.
# (purecn.example.output$candidates).
data(purecn.example.output)

# The max.candidate.solutions argument is set to a very low value only to 
# speed-up this example. This is not a good idea for real samples.
ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
    gatk.tumor.file=gatk.tumor.file, 
    vcf.file=vcf.file, genome="hg19",
    sampleid='Sample1', gc.gene.file=gc.gene.file, 
    candidates=purecn.example.output$candidates, max.candidate.solutions=2,
    fun.segmentation=segmentationCBS, args.segmentation=list(alpha=0.001))
}
