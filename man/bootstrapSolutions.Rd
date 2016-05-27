\name{bootstrapSolutions}
\alias{bootstrapSolutions}
\title{Filter unlikely purity/ploidy solutions}
\description{This function bootstraps SNVs, then re-ranks solutions 
by using the bootstrap estimate of the likelihood score, and then keeps 
only solutions that were ranked highest in any bootstrap replicate.
Large-scale copy number artifacts can cause true purity/ploidy 
solutions rank low.}
\usage{bootstrapSolutions(ret, n = 500, top = 2)}
\arguments{
  \item{ret}{Return object of the runAbsoluteCN() function.}
  \item{n}{Number of bootstrap replicates.}
  \item{top}{Include solution if it appears in the top n solutions of
any bootstrap replicate.}
}



\author{Markus Riester}




\examples{
data(purecn.example.output)
ret.boot <- bootstrapSolutions(purecn.example.output)
plotAbs(ret.boot, type="overview")
}
