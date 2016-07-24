\name{getSexFromCoverage}
\alias{getSexFromCoverage}
\title{Get sample sex from coverage}
\description{This function determines the sex of a sample by the coverage 
ratio of chrX and chrY. Loss of chromosome Y (LOY) can result in a wrong
female call. For small targeted panels, this will only work 
when sufficient sex marker genes such as AMELY are covered.
For optimal results, parameters might need to be tuned for the assay.}
\usage{getSexFromCoverage(gatk.coverage, min.ratio = 25, min.ratio.na = 20, 
    remove.outliers = TRUE, verbose = TRUE)}
\arguments{
  \item{gatk.coverage}{GATK coverage file or data read with \code{\link{readCoverageGatk}}.}
  \item{min.ratio}{Min chrX/chrY coverage ratio to call sample as female.}
  \item{min.ratio.na}{Min chrX/chrY coverage ratio to call sample as \code{NA}. 
This ratio defines a grey zone from \code{min.ratio.na} to 
\code{min.ratio} in which samples are not called.
The default is set to a copy number ratio that would be rare in male 
samples, but lower than expected in female samples. Contamination can be a 
source of ambiguous calls. Mappability issues on chromosome Y resulting in 
low coverage need to be considered when setting cutoffs.}
  \item{remove.outliers}{Removes coverage outliers before calculating mean chromosome coverages.}
  \item{verbose}{Verbose output.}
}

\value{Returns a \code{character(1)} with \code{M} for male, \code{F} 
for female, or \code{NA} if unknown.    }

\author{Markus Riester}



\seealso{\code{\link{getSexFromVcf}}}
\examples{
gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")
sex <- getSexFromCoverage(gatk.tumor.file)
}
