\name{autoCurateResults}
\alias{autoCurateResults}
\title{Heuristics to find the best purity/ploidy solution.}
\description{This implements a workflow with various heuristics, with the goal
of identifying correct purity/ploidy solutions in difficult samples.
This is mainly for automated copy number calling. This function 
may evolve over time and might produce different rankings
after PureCN updates.}
\usage{autoCurateResults(res, bootstrap = TRUE, verbose = TRUE)}
\arguments{
  \item{res}{Return object of the runAbsoluteCN() function.}
  \item{bootstrap}{Try to reduce the number of local optima by using the
bootstrapResults function.}
  \item{verbose}{Verbose output.}
}



\author{Markus Riester}





