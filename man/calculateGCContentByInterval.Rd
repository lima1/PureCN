\name{calculateGCContentByInterval}
\alias{calculateGCContentByInterval}
\title{Calculates GC content by interval}
\description{Uses \code{scanFa} from the Rsamtools package to retrieve GC 
content of intervals in a reference FASTA file.}
\usage{calculateGCContentByInterval(interval.file, reference.file, 
    output.file = NULL, ...)}
\arguments{
  \item{interval.file}{File specifying the intervals. Interval is expected in 
first column in format CHR:START-END. 
Instead of a file, a \code{GRanges} object can be provided.
This allows the use of BED files for example. Note that
GATK interval files are 1-based (first position of the genome
is 1). Other formats like BED files are often 0-based.
The \code{import} function will automatically convert to 
1-based \code{GRanges}.}
  \item{reference.file}{Reference FASTA file.}
  \item{output.file}{Optionally, write GC content file. }
  \item{\dots}{Additional parameters passed to the \code{read.delim}
function that reads the \code{interval.file}.}
}

\value{Returns GC content by interval.}

\author{Markus Riester}




\examples{
reference.file <- system.file("extdata", "ex2_reference.fa", 
    package="PureCN", mustWork = TRUE)
interval.file <- system.file("extdata", "ex2_intervals.txt", 
    package="PureCN", mustWork = TRUE)
bed.file <- system.file("extdata", "ex2_intervals.bed", 
    package="PureCN", mustWork = TRUE)
calculateGCContentByInterval(interval.file, reference.file, 
    output.file="gc_file.txt")

intervals <- import(bed.file)
calculateGCContentByInterval(intervals, reference.file, 
    output.file="gc_file.txt")
}
