\name{readCoverageGatk}
\alias{readCoverageGatk}
\title{readCoverageGatk}
\description{Read coverage file produced by The Genome Analysis Toolkit.
The three only important columns in the GATK-generated file
are: Target, total_coverage, and average_coverage. }
\usage{readCoverageGatk(file)}
\arguments{
  \item{file}{Exon coverage file as produced by GATK.}
}

\value{A data.frame with the parsed coverage information. }

\author{Markus Riester}




\examples{
gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")
coverage <- readCoverageGatk(gatk.tumor.file)
}
