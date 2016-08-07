\name{createExonWeightFile}
\alias{createExonWeightFile}
\title{Calculate exon weights}
\description{This function has been renamed to \code{\link{createTargetWeights}} 
and will be made defunct in the next Bioconductor release.}
\usage{createExonWeightFile(gatk.tumor.files, gatk.normal.files, 
    exon.weight.file)}
\arguments{
  \item{gatk.tumor.files}{A small number (1-3) of GATK tumor or normal coverage samples.}
  \item{gatk.normal.files}{A large number of GATK normal coverage samples (>20) 
to estimate exon log-ratio standard deviations.
Should not overlap with files in \code{gatk.tumor.files}.}
  \item{exon.weight.file}{Output filename.}
}



\author{Markus Riester}





