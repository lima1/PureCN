\name{createTargetWeights}
\alias{createTargetWeights}
\title{Calculate target weights}
\description{Creates a target weight file useful for segmentation. Requires a
set of GATK coverage files from normal samples. A small number of 
tumor (or other normal) samples is then tested against all normals. Target
weights will be set proportional to the inverse of coverage standard
deviation across all normals. Targets with high variance in coverage in the
pool of normals are thus down-weighted.}
\usage{createTargetWeights(gatk.tumor.files, gatk.normal.files, 
    target.weight.file, verbose = TRUE)}
\arguments{
  \item{gatk.tumor.files}{A small number (1-3) of GATK tumor or normal coverage samples.}
  \item{gatk.normal.files}{A large number of GATK normal coverage samples (>20) 
to estimate target log-ratio standard deviations.
Should not overlap with files in \code{gatk.tumor.files}.}
  \item{target.weight.file}{Output filename.}
  \item{verbose}{Verbose output.}
}

\value{A \code{data.frame} with target weights.}

\author{Markus Riester}




\examples{
target.weight.file <- "target_weights.txt"
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gatk.normal2.file <- system.file("extdata", "example_normal2.txt", 
    package="PureCN")
gatk.normal.files <- c(gatk.normal.file, gatk.normal2.file)
gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")

createTargetWeights(gatk.tumor.file, gatk.normal.files, target.weight.file)
}
