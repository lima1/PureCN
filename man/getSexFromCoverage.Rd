\name{getSexFromCoverage}
\alias{getSexFromCoverage}
\title{Get sample sex from coverage}
\description{This function determines the sex of a sample by the coverage 
ratio of chrX and chrY.}
\usage{getSexFromCoverage(gatk.coverage, min.ratio = 20, verbose = TRUE)}
\arguments{
  \item{gatk.coverage}{GATK coverage file or data read with readCoverageGatk.}
  \item{min.ratio}{Min chrX/chrY coverage ratio to call sample as female.}
  \item{verbose}{Verbose output.}
}

\value{Returns "M" for male, "F" for female, or NULL if unknown.    }

\author{Markus Riester}





