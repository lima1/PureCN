\name{createNormalDatabase}
\alias{createNormalDatabase}
\title{Create database of normal samples}
\description{Function to create a database of normal samples, used to find 
a good match for tumor copy number normalization.}
\usage{createNormalDatabase(gatk.normal.files, sex = NULL, 
    ...)}
\arguments{
  \item{gatk.normal.files}{Vector with file names pointing to GATK coverage files 
of normal samples. }
  \item{sex}{Vector of sex."F" for female, "M" for male. If NULL
determine from coverage.}
  \item{\dots}{Arguments passed to the prcomp function.}
}

\value{A normal database that can be used in the findBestNormal function to 
retrieve good matching normal samples for a given tumor sample.}

\author{Markus Riester}




\examples{
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gatk.normal2.file <- system.file("extdata", "example_normal2.txt", 
    package="PureCN")
gatk.normal.files <- c(gatk.normal.file, gatk.normal2.file)
normalDB <- createNormalDatabase(gatk.normal.files)
}
