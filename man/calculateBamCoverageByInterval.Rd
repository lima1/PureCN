\name{calculateBamCoverageByInterval}
\alias{calculateBamCoverageByInterval}
\title{Function to calculate coverage from BAM file}
\description{Takes a BAM file and an interval file as input and 
returns coverage for each interval. Coverage should be GC-normalized
using the \code{\link{correctCoverageBias}} function before determining
purity and ploidy with \code{\link{runAbsoluteCN}}.}
\usage{calculateBamCoverageByInterval(bam.file, interval.file, 
    output.file = NULL)}
\arguments{
  \item{bam.file}{Filename of a BAM file.}
  \item{interval.file}{File specifying the intervals. Interval is expected in 
first column in format CHR:START-END. The \code{gc.gene.file} can be used.}
  \item{output.file}{Optionally, write minimal coverage file. Can be read with the
\code{\link{readCoverageGatk}} function.}
}

\value{Returns total and average coverage by intervals.}

\author{Markus Riester}



\seealso{\code{\link{calculateGCContentByInterval} 
\link{correctCoverageBias} \link{runAbsoluteCN}}}
\examples{
bam.file <- system.file("extdata", "ex1.bam", package="PureCN", 
    mustWork = TRUE)
interval.file <- system.file("extdata", "ex1_intervals.txt", 
    package="PureCN", mustWork = TRUE)

coverage <- calculateBamCoverageByInterval(bam.file=bam.file, 
    interval.file=interval.file)
}
