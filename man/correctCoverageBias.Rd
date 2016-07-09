\name{correctCoverageBias}
\alias{correctCoverageBias}
\title{Correct for GC bias}
\description{Takes as input coverage data in GATK format (or data 
read by \code{\link{readCoverageGatk}}) and a mapping file for 
GC content, and then uses a loess normalization for bias correction. 
Largely follows the GC correction of the TitanCNA package.}
\usage{correctCoverageBias(gatk.coverage.file, gc.gene.file, 
    output.file = NULL)}
\arguments{
  \item{gatk.coverage.file}{Exon coverage file as produced by GATK. Either a file name
or data parsed with the readCoverageGatk function.}
  \item{gc.gene.file}{File providing GC content for each exon in the coverage files.
First column in format CHR:START-END. Second column GC content (0 to 1). 
Third column provides gene symbols, which are optional, but used in
\code{\link{runAbsoluteCN}} to generate gene level calls. This file
can be generated with GATK GCContentByInterval tool or with the
\code{\link{calculateGCContentByInterval}} function.}
  \item{output.file}{Optionally, write file with GC corrected coverage. Can be read with
the \code{\link{readCoverageGatk}} function.}
}

\value{GC normalized coverage.}

\author{Markus Riester}



\seealso{\code{\link{calculateGCContentByInterval}}}
\examples{
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
    package="PureCN")
coverage <- correctCoverageBias(gatk.normal.file, gc.gene.file)
}
