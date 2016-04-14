\name{chr.hash}
\alias{chr.hash}
\docType{data}
\title{
A data.frame of chromosome names to hash to numeric representation.
}
\description{
a table of chromosome name and their corresponding numerical representation. This is needed because the way DNAcopy encodes the chromosomes.
}
\usage{data(chr.hash)}
\value{
    A data frame with 24 observations on the following 2 variables.
    \describe{
        \item{\code{chr}}{a factor with levels \code{chr1} \code{chr10} \code{chr11} \code{chr12} \code{chr13} \code{chr14} \code{chr15} \code{chr16} \code{chr17} \code{chr18} \code{chr19} \code{chr2} \code{chr20} \code{chr21} \code{chr22} \code{chr3} \code{chr4} \code{chr5} \code{chr6} \code{chr7} \code{chr8} \code{chr9} \code{chrX} \code{chrY}}
        \item{\code{number}}{a numeric vector}
    }
}
\references{
https://secure.genome.ucla.edu/index.php/ExomeCNV_User_Guide
}
\examples{
data(chr.hash)
}
\keyword{datasets}
