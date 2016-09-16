\name{calculateLogRatio}
\alias{calculateLogRatio}
\title{Calculate coverage log-ratio of tumor vs. normal}
\description{This function is automatically called by \code{\link{runAbsoluteCN}} 
when normal and tumor coverage are provided (and not a segmentation file 
or target-level log-ratios). This function is therefore 
normally not called by the user.}
\usage{calculateLogRatio(normal, tumor, verbose = TRUE)}
\arguments{
  \item{normal}{Normal coverage read in by the 
\code{\link{readCoverageGatk}} function.}
  \item{tumor}{Tumor coverage read in by the 
\code{\link{readCoverageGatk}} function.}
  \item{verbose}{Verbose output.}
}



\author{Markus Riester}




\examples{
normal.coverage.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")
normal <- readCoverageGatk(normal.coverage.file)
tumor <- readCoverageGatk(tumor.coverage.file)
log.ratio <- calculateLogRatio(normal, tumor, verbose=FALSE)
}
