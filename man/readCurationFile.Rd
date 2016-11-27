\name{readCurationFile}
\alias{readCurationFile}
\title{Read curation file}
\description{Function that can be used to read the curated output of
the \code{\link{runAbsoluteCN}} function.}
\usage{readCurationFile(file.rds, file.curation = gsub(".rds$", 
    ".csv", file.rds), remove.failed = FALSE, report.best.only = FALSE, 
    min.ploidy = NULL, max.ploidy = NULL, verbose = FALSE)}
\arguments{
  \item{file.rds}{Output of the \code{\link{runAbsoluteCN}} function, serialized with 
\code{saveRDS}.}
  \item{file.curation}{Filename of a curation file that points to the correct 
tumor purity and ploidy solution.}
  \item{remove.failed}{Do not return solutions that failed.}
  \item{report.best.only}{Only return correct/best solution (useful on low memory
machines when lots of samples are loaded).}
  \item{min.ploidy}{Minimum ploidy to be considered. If \code{NULL}, all. Can be 
used to automatically ignore unlikely solutions.}
  \item{max.ploidy}{Maximum ploidy to be considered. If \code{NULL}, all. Can be 
used to automatically ignore unlikely solutions.}
  \item{verbose}{Verbose output.}
}

\value{The return value of the corresponding \code{\link{runAbsoluteCN}} call,
but with the results array manipulated according the curation CSV file 
and arguments of this function.}

\author{Markus Riester}



\seealso{\code{\link{runAbsoluteCN} \link{createCurationFile}}}
\examples{
data(purecn.example.output)
file.rds <- 'Sample1_PureCN.rds'
createCurationFile(file.rds) 
# User can change the maximum likelihood solution manually in the generated 
# CSV file. The correct solution is then loaded with readCurationFile.
purecn.curated.example.output <-readCurationFile(file.rds) 
}
