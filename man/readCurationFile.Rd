\name{readCurationFile}
\alias{readCurationFile}
\title{readCurationFile}
\description{Function that can be used to read the curated output of the runAbsoluteCN function.}
\usage{readCurationFile(file.rds, file.curation = gsub(".rds$", 
    ".csv", file.rds), remove.failed = FALSE, report.best.only = FALSE, 
    min.ploidy = NULL, max.ploidy = NULL)}
\arguments{
  \item{file.rds}{Output of the runAbsoluteCN() function, serialized with saveRDS()}
  \item{file.curation}{Filename that points to the correct solution}
  \item{remove.failed}{Do not return solutions that failed.}
  \item{report.best.only}{Only return correct/best solution (useful on low memory machines when lots of samples are loaded)}
  \item{min.ploidy}{Minimum ploidy to be considered. If NULL, all. Can be used to automatically ignore unlikely solutions.}
  \item{max.ploidy}{Maximum ploidy to be considered. If NULL, all. Can be used to automatically ignore unlikely solutions.}
}

\value{The return value of the corresponding runAbsoluteCN call, but with the
results array manipulated according the curation CSV file and arguments
of this function.}

\author{Markus Riester}




\examples{
data(purecn.example.output)
file.rds <- 'Sample1_PureCN.rds'
createCurationFile(file.rds) 
# User can change the maximum likelihood solution manually in the generated 
# CSV file. The correct solution is then loaded with readCurationFile.
purecn.curated.example.output <-readCurationFile(file.rds) 
}
