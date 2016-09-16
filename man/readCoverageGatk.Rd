\name{readCoverageGatk}
\alias{readCoverageGatk}
\title{Read GATK coverage files}
\description{Read coverage file produced by The Genome Analysis Toolkit or
by \code{\link{calculateBamCoverageByInterval}}.}
\usage{readCoverageGatk(file)}
\arguments{
  \item{file}{Exon coverage file as produced by GATK.}
}

\value{A \code{data.frame} with the parsed coverage information. }

\author{Markus Riester}



\seealso{\code{\link{calculateBamCoverageByInterval}}}
\examples{
tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")
coverage <- readCoverageGatk(tumor.coverage.file)
}
