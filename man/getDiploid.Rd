\name{getDiploid}
\alias{getDiploid}
\title{Function to extract diploid solutions.}
\description{This function can be used to extract purity and ploidy solutions that
are diploid with only few CNVs. Since high ploidy solutions 
typically have a very different copy number profile, one of the 
identified diploid solutions is likely correct if there are any. This 
function can be used for automated curation workflows; very silent genomes
have by definition only a small number of CNVs, making it difficult for the
algorithm to correctly identify purity and ploidy. If the maximum
likelihood solution is diploid, it is always returned; all other solutions
must pass the more stringent criteria as defined in the function arguments.}
\usage{getDiploid(res, min.diploid = 0.5, min.single.gain.loss = 0.05, 
    max.non.single.gain.loss = 0.1, max.loh = 0.5, 
    min.log.likelihood = NULL)}
\arguments{
  \item{res}{Return object of the \code{\link{runAbsoluteCN}} function.}
  \item{min.diploid}{Minimum fraction of genome with normal copy number 2.}
  \item{min.single.gain.loss}{Minimum fraction of genome with copy number 1 or 3. This makes sure that 
low purity samples are not confused with quiet samples.}
  \item{max.non.single.gain.loss}{Maximum fraction of genome with copy number smaller 1 or more than 3.}
  \item{max.loh}{Maximum fraction of genome in LOH. }
  \item{min.log.likelihood}{Minimum copy number log-likelihood to consider sample. If 
\code{NULL}, not tested.}
}

\value{A list with elements
\item{ids}{The ids of diploid solutions
(\code{res$results[ids]}).}
\item{fraction.non.single}{The fraction of the genome
with copy number <1 or >3.}}

\author{Markus Riester}



\seealso{\code{\link{runAbsoluteCN}}}
\examples{
data(purecn.example.output)
# no diploid solutions in the example
idx <- getDiploid(purecn.example.output)
}
