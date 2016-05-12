\name{calculateBamCoverageByInterval}
\alias{calculateBamCoverageByInterval}
\title{Function to calculate coverage from BAM file}
\description{Takes a BAM file and an interval file as input and 
returns coverage for each interval.}
\usage{calculateBamCoverageByInterval(bam.file, interval.file, 
    output.file = NULL)}
\arguments{
  \item{bam.file}{Filename of a BAM file.}
  \item{interval.file}{File specifying the interval. Interval is expected in 
first column in format CHR:START-END. The gc.gene.file can be used.}
  \item{output.file}{Optionally, write minimal coverage file. Can be read with the
readCoverageGatk function.}
}

\value{Returns total and average coverage by intervals.}

\author{Markus Riester}




\examples{
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
    package="PureCN")
# calculateBamCoverageByInterval(bam.file, gc.gene.file)   
}
