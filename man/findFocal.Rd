\name{findFocal}
\alias{findFocal}
\title{Find focal amplifications}
\description{Function to find focal amplifications in segmented data. 
This is automatically called in runAbsoluteCN.}
\usage{findFocal(seg, size.cutoff = 2e+06, cn.diff = 2, amp.cutoff = 6)}
\arguments{
  \item{seg}{Segmentation data.}
  \item{size.cutoff}{Cutoff for focal in base pairs.}
  \item{cn.diff}{Minimum copy number delta between neighboring segments.}
  \item{amp.cutoff}{Minimum amplification integer copy number.}
}

\value{Boolean vector for each segment wether it is focally 
amplified or not.}

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

# Speed-up the runAbsoluteCN call by using the stored grid-search 
# (purecn.example.output$candidates).
data(purecn.example.output)

# The max.candidate.solutions parameter is set to a very low value only to
# speed-up this example.  This is not a good idea for real samples.
ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
    gatk.tumor.file=gatk.tumor.file,
    vcf.file=vcf.file, sampleid='Sample1', gc.gene.file=gc.gene.file,
    candidates=purecn.example.output$candidates, max.candidate.solutions=2,
    args.focal=list(size.cutoff = 2e+06), fun.focal=findFocal)
}
