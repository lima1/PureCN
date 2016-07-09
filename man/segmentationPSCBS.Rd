\name{segmentationPSCBS}
\alias{segmentationPSCBS}
\title{PSCBS segmentation}
\description{Segmentation function. Uses the PSCBS package.
This function is called via the fun.segmentation argument of runAbsoluteCN.
The arguments are passed via args.segmentation.}
\usage{segmentationPSCBS(normal, tumor, log.ratio, plot.cnv, 
    coverage.cutoff, sampleid = sampleid, exon.weight.file = NULL, 
    flavor = "tcn&dh", tauA = 0.03, vcf = NULL, tumor.id.in.vcf = 1, 
    normal.id.in.vcf = NULL, max.segments = NULL, chr.hash = NULL, 
    verbose = TRUE, ...)}
\arguments{
  \item{normal}{GATK coverage data for normal sample.}
  \item{tumor}{GATK coverage data for tumor sample.}
  \item{log.ratio}{Copy number log-ratios, one for each exon in coverage file.}
  \item{plot.cnv}{Segmentation plots.}
  \item{coverage.cutoff}{Minimum coverage in both normal and tumor.}
  \item{sampleid}{Sample id, used in output files.}
  \item{exon.weight.file}{Can be used to assign weights to exons. NOT SUPPORTED YET.}
  \item{flavor}{Flavor value for PSCBS. See segmentByNonPairedPSCBS.}
  \item{tauA}{tauA argument for PSCBS. See segmentByNonPairedPSCBS.}
  \item{vcf}{Optional VCF object with germline allelic ratios.}
  \item{tumor.id.in.vcf}{Id of tumor in case multiple samples are stored in VCF.}
  \item{normal.id.in.vcf}{Id of normal in in VCF. If NULL, use unpaired PSCBS.}
  \item{max.segments}{If not NULL, try a higher undo.SD parameter if number of
segments exceeds the threshold.}
  \item{chr.hash}{Mapping of non-numerical chromsome names to numerical names
(e.g. chr1 to 1, chr2 to 2, etc.). If NULL, assume chromsomes
are properly ordered.   }
  \item{verbose}{Verbose output.}
  \item{\dots}{Additional parameters passed to the segmentByNonPairedPSCBS function.}
}

\value{A list with elements seg and size. "seg" contains the 
segmentation, "size" the size of all segments in base pairs.    }

\author{Markus Riester}




\examples{
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")
vcf.file <- system.file("extdata", "example_vcf.vcf", 
    package="PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
    package="PureCN")

# The max.candidate.solutions argument is set to a very low value only to 
# speed-up this example. This is not a good idea for real samples.
#ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
#    gatk.tumor.file=gatk.tumor.file, vcf.file=vcf.file, 
#    genome="hg19", sampleid='Sample1', 
#    max.candidate.solutions=2, 
#    gc.gene.file=gc.gene.file, fun.segmentation=segmentationPSCBS)
}
