\name{filterTargets}
\alias{filterTargets}
\title{Remove low quality targets}
\description{This function determines which intervals in the coverage files should
be included or excluded in the segmentation. It is called via the
\code{fun.filterTargets} argument of \code{\link{runAbsoluteCN}}. The 
arguments are passed via \code{args.filterTargets}.}
\usage{filterTargets(log.ratio, tumor, gc.data, seg.file, 
    filter.lowhigh.gc = 0.001, filter.targeted.base = 4, 
    normalDB = NULL, normalDB.min.coverage = 0.2, verbose)}
\arguments{
  \item{log.ratio}{Copy number log-ratios, one for each target or interval
in coverage file.}
  \item{tumor}{GATK coverage data for tumor sample.}
  \item{gc.data}{\code{data.frame} with GC bias for each interval.}
  \item{seg.file}{If not \code{NULL}, then do not filter targets, because data is
already segmented via the provided segmentation file.}
  \item{filter.lowhigh.gc}{Quantile q (defines lower q and upper 1-q) 
for removing targets with outlier GC profile. Assuming that GC correction 
might not have been worked on those. Requires \code{gc.gene.file}.}
  \item{filter.targeted.base}{Exclude exons with targeted base (size) smaller 
than this cutoff. This is useful when the same interval file was used to
calculate GC content. For such small exons, the GC content is likely 
very different from the true GC content of the probes.}
  \item{normalDB}{Normal database, created with \code{\link{createNormalDatabase}}.}
  \item{normalDB.min.coverage}{Exclude targets with coverage lower than 20 percent of the 
chromosome median in the pool of normals.}
  \item{verbose}{Verbose output.}
}



\author{Markus Riester}




\examples{
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gatk.normal2.file <- system.file("extdata", "example_normal2.txt", 
    package="PureCN")
gatk.normal.files <- c(gatk.normal.file, gatk.normal2.file)
normalDB <- createNormalDatabase(gatk.normal.files)

gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")
vcf.file <- system.file("extdata", "example_vcf.vcf", 
    package="PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
    package="PureCN")

# Speed-up the runAbsoluteCN call by using the stored grid-search 
# (purecn.example.output$candidates).
data(purecn.example.output)

# The max.candidate.solutions, candidates, max.ploidy and test.purity 
# parameters are set to non-default values to speed-up this example.  
# This is not a good idea for real samples.

ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
    gatk.tumor.file=gatk.tumor.file, 
    genome="hg19", vcf.file=vcf.file, sampleid='Sample1', 
    gc.gene.file=gc.gene.file,
    args.filterTargets=list(normalDB=normalDB),
    candidates=purecn.example.output$candidates, 
    max.ploidy=4, test.purity=seq(0.3, 0.7, by=0.05),
    max.candidate.solutions=1)
}
