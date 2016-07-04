\name{centromeres}
\alias{centromeres}
\docType{data}
\title{
A list of data.frames containing centromere positions.
}
\description{
A list of data.frames containing centromere positions for hg18, hg19 and hg38.
Downloaded from the UCSC genome browser.
}
\usage{data(centromeres)}
\value{
    A list with three data frames, "hg18", "hg19", and "hg38". Each containes
three columns
    \describe{
        \item{\code{chrom}}{a factor with levels \code{chr1} \code{chr10} \code{chr11} \code{chr12} \code{chr13} \code{chr14} \code{chr15} \code{chr16} \code{chr17} \code{chr18} \code{chr19} \code{chr2} \code{chr20} \code{chr21} \code{chr22} \code{chr3} \code{chr4} \code{chr5} \code{chr6} \code{chr7} \code{chr8} \code{chr9} \code{chrX} \code{chrY}}
        \item{\code{chromStart}}{a numeric vector}
        \item{\code{chromEnd}}{a numeric vector}
    }
}
\references{
The script downloadCentromeres.R in the extdata directory was used to generate
the data.frames.
}
\examples{
data(centromeres)
}
\keyword{datasets}
