\name{calculateGCContentByInterval}
\alias{calculateGCContentByInterval}
\title{Calculates GC content by interval}
\description{Uses \code{scanFa} from the Rsamtools package to retrieve GC 
content of intervals in a reference FASTA file.}
\usage{calculateGCContentByInterval(interval.file, reference.file, 
    output.file = NULL)}
\arguments{
  \item{interval.file}{File specifying the intervals. Interval is expected in 
first column in format CHR:START-END.}
  \item{reference.file}{Reference FASTA file.}
  \item{output.file}{Optionally, write GC content file. }
}

\value{Returns GC content by interval.}

\author{Markus Riester}




\examples{
reference.file <- system.file("extdata", "ex2_reference.fa", 
    package="PureCN", mustWork = TRUE)
interval.file <- system.file("extdata", "ex2_intervals.txt", 
    package="PureCN", mustWork = TRUE)
calculateGCContentByInterval(interval.file, reference.file, 
    output.file="gc_file.txt")
}
