\name{createCurationFile}
\alias{createCurationFile}
\title{Create file to curate PureCN results}
\description{Function to create a CSV file that can be used to mark the correct solution
in the output of a runAbsoluteCN() run.}
\usage{createCurationFile(file.rds, overwrite.uncurated = TRUE)}
\arguments{
  \item{file.rds}{Output of the runAbsoluteCN() function, serialized with saveRDS()}
  \item{overwrite.uncurated}{Overwrite existing files unless flagged as "Curated".}
}

\value{A data.frame with the tumor purity and ploidy of the maximum likelihood
solution }

\author{Markus Riester}




\examples{
data(purecn.example.output)
file.rds <- 'Sample1_PureCN.rds'
saveRDS(purecn.example.output, file=file.rds)
createCurationFile(file.rds) 
}
