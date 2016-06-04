\name{segmentationPSCBS}
\alias{segmentationPSCBS}
\title{PSCBS segmentation}
\description{Segmentation function. Uses the PSCBS package.
This function is called via the fun.segmentation argument of runAbsoluteCN.
The arguments are passed via args.segmentation.}
\usage{segmentationPSCBS(normal, tumor, log.ratio, plot.cnv, 
    coverage.cutoff, sampleid = sampleid, exon.weight.file = NULL, 
    flavor = "tcn&dh", tauA = 0.03, vcf = NULL, tumor.id.in.vcf = 1, 
    verbose = TRUE, ...)}
\arguments{
  \item{normal}{GATK coverage data for normal sample.}
  \item{tumor}{GATK coverage data for tumor sample.}
  \item{log.ratio}{Copy number log-ratios, one for each exon in coverage file.}
  \item{plot.cnv}{Segmentation plots.}
  \item{coverage.cutoff}{Minimum coverage in both normal and tumor.}
  \item{sampleid}{Sample id, used in output files.}
  \item{exon.weight.file}{Can be used to assign weights to exons. NOT SUPPORTED YET.}
  \item{flavor}{Flavor value for PSBCS. See segmentByNonPairedPSCBS.}
  \item{tauA}{tauA argument for PSCBS. See segmentByNonPairedPSCBS.}
  \item{vcf}{Optional VCF object with germline allelic ratios.}
  \item{tumor.id.in.vcf}{Id of tumor in case multiple samples are stored in VCF.}
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
ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
    gatk.tumor.file=gatk.tumor.file, vcf.file=vcf.file, sampleid='Sample1', 
    max.candidate.solutions=2, remove.off.target.snvs=TRUE,
    gc.gene.file=gc.gene.file, fun.segmentation=segmentationPSCBS)
}
