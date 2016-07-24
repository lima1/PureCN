\name{poolCoverage}
\alias{poolCoverage}
\title{Pool coverage from multiple samples}
\description{Averages the coverage of a list of samples.}
\usage{poolCoverage(all.data, remove.chrs = c(), w = NULL)}
\arguments{
  \item{all.data}{List of normals, read with \code{\link{readCoverageGatk}}.}
  \item{remove.chrs}{Remove these chromosomes from the pool.}
  \item{w}{List of weights for each sample.}
}

\value{A \code{data.frame} with the averaged coverage over all normals.}

\author{Markus Riester}



\seealso{\code{\link{readCoverageGatk}}}
\examples{
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gatk.normal2.file <- system.file("extdata", "example_normal2.txt", 
    package="PureCN")
gatk.normal.files <- c(gatk.normal.file, gatk.normal2.file)
gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")

normalDB <- createNormalDatabase(gatk.normal.files)

# get the best 2 normals and average them
gatk.best.normal.files <- findBestNormal(gatk.tumor.file, normalDB, 
    num.normals=2)
pool <- poolCoverage(lapply(gatk.best.normal.files, readCoverageGatk),
     remove.chrs=c('chrX', 'chrY'))
}
