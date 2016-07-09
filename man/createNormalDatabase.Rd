\name{createNormalDatabase}
\alias{createNormalDatabase}
\title{Create database of normal samples}
\description{Function to create a database of normal samples, used to find 
a good match for tumor copy number normalization. Internally, this 
function determines the sex of the samples and trains a PCA
that is later used for clustering a tumor file with all normal samples 
in the database.}
\usage{createNormalDatabase(gatk.normal.files, sex = NULL, 
    ...)}
\arguments{
  \item{gatk.normal.files}{Vector with file names pointing to GATK coverage files 
of normal samples. }
  \item{sex}{Vector of sex."F" for female, "M" for male. If all chromosomes are diploid, specify "diploid". 
If NULL determine from coverage.}
  \item{\dots}{Arguments passed to the \code{prcomp} function.}
}

\value{A normal database that can be used in the \code{\link{findBestNormal}} 
function to retrieve good matching normal samples for a given tumor sample.}

\author{Markus Riester}



\seealso{\code{\link{findBestNormal}}}
\examples{
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gatk.normal2.file <- system.file("extdata", "example_normal2.txt", 
    package="PureCN")
gatk.normal.files <- c(gatk.normal.file, gatk.normal2.file)
normalDB <- createNormalDatabase(gatk.normal.files)
}
