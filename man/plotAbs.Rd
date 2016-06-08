\name{plotAbs}
\alias{plotAbs}
\title{Plots for analyzing PureCN solutions}
\description{This function provides various plots for finding correct 
purity and ploidy combinations in the results of a runAbsoluteCN call.}
\usage{plotAbs(res, ids = NULL, type = c("hist", "overview", 
    "BAF", "AF", "LOH", "all"), chr = NULL, germline.only = TRUE, 
    show.contour = FALSE, purity = NULL, ploidy = NULL, 
    alpha = TRUE, show.segment.means = c("SNV", "segments", 
        "both"), ...)}
\arguments{
  \item{res}{Return object of the runAbsoluteCN() function.}
  \item{ids}{Candidate solutions to be plotted. ids=1 will draw the 
plot for the maximum likelihood solution.}
  \item{type}{Different types of plots. "hist" will plot a histogram, 
assigning log-ratio peaks to integer values. "overview" will plot all 
local optima, sorted by likelihood. "BAF" plots something like a B-allele
frequency plot known from SNP arrays: it plots allele frequencies of 
germline variants (or most likely germline when status is not available) 
against copy number. AF plots observed allelic fractions against expected
(purity), maximum likelihood (optimal multiplicity) allelic fractions. 
"all" plots all, and is useful for generate a PDF for a sample for manual
inspection.}
  \item{chr}{If NULL, show all chromosomes, otherwise only the ones 
specified (type=BAF only).}
  \item{germline.only}{If TRUE, show only variants most likely being germline in 
BAF plot. Useful to set to FALSE (in combination with chr) to study 
potential artifacts.}
  \item{show.contour}{For type overview, display contour plot.}
  \item{purity}{Display expected integer copy numbers for purity, defaults 
to purity of the solution (type=hist only).}
  \item{ploidy}{Display expected integer copy numbers for ploidy, defaults 
to ploidy of the solution (type=hist only).}
  \item{alpha}{Add transparency to the plot if VCF contains many variants 
(>2000, type=AF and BAF only). }
  \item{show.segment.means}{Show segment means in germline allele frequency plot? 
If "both", show SNVs and segment means. If "SNV" show
all SNVs. Type "AF" only.}
  \item{\dots}{Additonal parameters passed to the plot() function. }
}

\value{Returns NULL}

\author{Markus Riester}




\examples{
data(purecn.example.output)
plotAbs(purecn.example.output, type="overview")
# plot details for the maximum likelihood solution (rank 1)
plotAbs(purecn.example.output, 1, type="hist")
plotAbs(purecn.example.output, 1, type="BAF")
}
