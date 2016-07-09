\name{createExonWeightFile}
\alias{createExonWeightFile}
\title{Calculate exon weights}
\description{Creates an exon weight file useful for segmentation. Requires a
set of GATK coverage files from normal samples. A small number of 
tumor (or other normal) samples is then tested against all normals. Exon
weights will be set proportional to the inverse of coverage standard
deviation across all normals. Exons with high variance in coverage in the
pool of normals are thus down-weighted.}
\usage{createExonWeightFile(gatk.tumor.files, gatk.normal.files, 
    exon.weight.file)}
\arguments{
  \item{gatk.tumor.files}{A small number (1-3) of GATK tumor or normal coverage samples.}
  \item{gatk.normal.files}{A large number of GATK normal coverage samples (>20) 
to estimate exon log-ratio standard deviations.
Should not overlap with files in \code{gatk.tumor.files}.}
  \item{exon.weight.file}{Output filename.}
}

\value{A data.frame with exon weights.}

\author{Markus Riester}




\examples{
exon.weight.file <- "exon_weights.txt"
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gatk.normal2.file <- system.file("extdata", "example_normal2.txt", 
    package="PureCN")
gatk.normal.files <- c(gatk.normal.file, gatk.normal2.file)
gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")

createExonWeightFile(gatk.tumor.file, gatk.normal.files, exon.weight.file)
}
