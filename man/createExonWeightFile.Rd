\name{createExonWeightFile}
\alias{createExonWeightFile}
\title{Calculate exon weights}
\description{Creates an exon weight file useful for segmentation, by down-weighting unreliable exons.}
\usage{createExonWeightFile(gatk.tumor.files, gatk.normal.files, 
    exon.weight.file)}
\arguments{
  \item{gatk.tumor.files}{A small number (1-3) of GATK tumor coverage samples.}
  \item{gatk.normal.files}{A large number of GATK normal coverage samples (>20) to estimate exon log-ratio standard deviations.}
  \item{exon.weight.file}{Output filename.}
}

\value{A data.frame with exon weights.}

\author{Markus Riester}




\examples{
exon.weight.file <- "exon_weights.txt"
gatk.normal.file <- system.file("extdata", "example_normal.txt", package="PureCN")
gatk.normal2.file <- system.file("extdata", "example_normal2.txt", package="PureCN")
gatk.normal.files <- c(gatk.normal.file, gatk.normal2.file)
gatk.tumor.file <- system.file("extdata", "example_tumor.txt", package="PureCN")

createExonWeightFile(gatk.tumor.file, gatk.normal.files, exon.weight.file)
}
