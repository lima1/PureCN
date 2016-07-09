\name{predictSomatic}
\alias{predictSomatic}
\title{Predict germline vs. somatic status}
\description{This function takes as input the output of a 
\code{\link{runAbsoluteCN}} run and provides SNV posterior probabilities
for all possible states. }
\usage{predictSomatic(res, id = 1, cutoff = 0.1)}
\arguments{
  \item{res}{Return object of the \code{\link{runAbsoluteCN}} function.}
  \item{id}{Candidate solutions to be analyzed. \code{id=1} will analyze the 
maximum likelihood solution.}
  \item{cutoff}{Exclude maternal/paternal chromosome number in segment if 
posterior probability is lower than cutoff.}
}

\value{A data.frame with SNV state posterior probabilities.    }

\author{Markus Riester}



\seealso{\code{\link{runAbsoluteCN}}}
\examples{
data(purecn.example.output)
# the output data was created using a matched normal sample, but in case
# no matched normal is available, this will help predicting somatic vs. 
# germline status
purecn.snvs <- predictSomatic(purecn.example.output)
}
