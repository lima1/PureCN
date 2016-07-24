\name{bootstrapResults}
\alias{bootstrapResults}
\title{Filter unlikely purity/ploidy solutions}
\description{This function bootstraps SNVs, then re-ranks solutions 
by using the bootstrap estimate of the likelihood score, and then keeps 
only solutions that were ranked highest in any bootstrap replicate.
Large-scale copy number artifacts can cause true purity/ploidy 
solutions rank low.}
\usage{bootstrapResults(res, n = 500, top = 2)}
\arguments{
  \item{res}{Return object of the \code{\link{runAbsoluteCN}} function.}
  \item{n}{Number of bootstrap replicates.}
  \item{top}{Include solution if it appears in the top \code{n} solutions of
any bootstrap replicate.}
}

\value{Returns the \code{\link{runAbsoluteCN}} object with low 
likelihood solutions removed. Also adds a bootstrap value 
to each solution. This value is the fraction of bootstrap replicates 
in which the solution ranked first.                }

\author{Markus Riester}



\seealso{\code{\link{runAbsoluteCN}}}
\examples{
data(purecn.example.output)
ret.boot <- bootstrapResults(purecn.example.output, n=100)
plotAbs(ret.boot, type="overview")
}
