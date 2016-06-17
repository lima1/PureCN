\name{callAlterations}
\alias{callAlterations}
\title{Calling of amplifications and deletions}
\description{Function to extract major copy number alterations from a 
runAbsoluteCN return object.}
\usage{callAlterations(res, id = 1, cutoffs = c(0.5, 6, 7), 
    log.ratio.cutoffs = c(-0.9, 0.9), failed = NULL, 
    ...)}
\arguments{
  \item{res}{Return object of the runAbsoluteCN() function.}
  \item{id}{Candidate solutions to be used. id=1 will use the 
maximum likelihood (or curated) solution.}
  \item{cutoffs}{Copy numbers cutoffs to call losses, focal amplifications 
and broad amplifications.}
  \item{log.ratio.cutoffs}{Copy numbers log-ratio cutoffs to call losses and amplifications 
in failed samples.}
  \item{failed}{Indicates whether sample was failed. If NULL, use available 
annotation, which can be set in the curation file.}
  \item{\dots}{
}
}

\value{A data.frame with gene-level amplification and deletion calls.}

\author{Markus Riester}




\examples{
data(purecn.example.output)
callAlterations(purecn.example.output)
}
