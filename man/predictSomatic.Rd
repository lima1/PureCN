\name{predictSomatic}
\alias{predictSomatic}
\title{Predict germline vs. somatic status}
\description{This function takes as input the output of a 
\code{\link{runAbsoluteCN}} run and provides SNV posterior probabilities
for all possible states. }
\usage{predictSomatic(res, id = 1, return.vcf = FALSE)}
\arguments{
  \item{res}{Return object of the \code{\link{runAbsoluteCN}} function.}
  \item{id}{Candidate solutions to be analyzed. \code{id=1} will analyze the 
maximum likelihood solution.}
  \item{return.vcf}{Returns an annotated \code{CollapsedVCF} object. Note that 
this VCF will only contain variants not filtered out by the 
\code{filterVcf} functions. Variants outside segments or intervals
might be included or not depending on \code{\link{runAbsoluteCN}}
arguments.}
}

\value{A \code{data.frame} or \code{CollapsedVCF} with SNV state 
posterior probabilities.    }

\author{Markus Riester}



\seealso{\code{\link{runAbsoluteCN}}}
\examples{
data(purecn.example.output)
# the output data was created using a matched normal sample, but in case
# no matched normal is available, this will help predicting somatic vs. 
# germline status
purecn.snvs <- predictSomatic(purecn.example.output)

# write a VCF file
purecn.vcf <- predictSomatic(purecn.example.output, return.vcf=TRUE)
writeVcf(purecn.vcf, file="Sample1_PureCN.vcf")
}
