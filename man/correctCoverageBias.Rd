\name{correctCoverageBias}
\alias{correctCoverageBias}
\title{Correct for GC bias}
\description{Takes as input coverage data in GATK format (or data read by readCoverageGatk)
and a mapping file for GC content, and then uses a loess normalization for
bias correction. Largely follows the GC correction of the TitanCNA package.}
\usage{correctCoverageBias(gatk.coverage.file, gc.gene.file, 
    output.file = NULL)}
\arguments{
  \item{gatk.coverage.file}{Exon coverage file as produced by GATK.}
  \item{gc.gene.file}{GC gene file, which provides GC content for each exon in the coverage files. 
First column in format CHR:START-END. Second column GC content (0 to 1). 
Third column provides gene symbols, which are optional, but used in
runAbsoluteCN to generate gene level calls.}
  \item{output.file}{Optionally, write minimal GATK coverage file with GC corrected coverage.}
}

\value{GC normalized coverage}

\author{Markus Riester}




\examples{
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
    package="PureCN")
coverage <- correctCoverageBias(gatk.normal.file, gc.gene.file)
}
