\name{findBestNormal}
\alias{findBestNormal}
\title{findBestNormal}
\description{Function to find the best matching normal for a provided tumor sample.}
\usage{findBestNormal(gatk.tumor.file, normalDB, pcs = 1:3, 
    num.normals = 1, ignore.sex = FALSE, verbose = TRUE)}
\arguments{
  \item{gatk.tumor.file}{GATK coverage file of a tumor sample.}
  \item{normalDB}{Database of normal samples, created with createNormalDatabase().}
  \item{pcs}{Principal components to use for distance calculation.}
  \item{num.normals}{Return the num.normals best normals.}
  \item{ignore.sex}{If FALSE, detects sex of sample returns a best normal 
with matching sex.}
  \item{verbose}{Verbose output.}
}

\value{Return the filename of the best matching normal    }

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
gatk.best.normal.file <- findBestNormal(gatk.tumor.file, normalDB)
}
