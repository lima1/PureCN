\name{filterVcfMuTect}
\alias{filterVcfMuTect}
\title{Filter VCF MuTect}
\description{Function to remove artifacts and low confidence/quality calls from
a MuTect generated VCF file. Also applies filters defined in
\code{filterVcfBasic}. This function will only keep variants 
listed in the stats file and those not matching the specified
failure reasons. }
\usage{filterVcfMuTect(vcf, tumor.id.in.vcf = NULL, stats.file = NULL, 
    ignore = c("clustered_read_position", "fstar_tumor_lod", 
        "nearby_gap_events", "poor_mapping_region_alternate_allele_mapq", 
        "poor_mapping_region_mapq0", "possible_contamination", 
        "strand_artifact", "seen_in_panel_of_normals"), 
    verbose = TRUE, ...)}
\arguments{
  \item{vcf}{\code{CollapsedVCF} object, read in with the \code{readVcf} function 
from the VariantAnnotation package.}
  \item{tumor.id.in.vcf}{The tumor id in the VCF file, optional.}
  \item{stats.file}{MuTect stats file}
  \item{ignore}{
}
  \item{verbose}{Verbose output.}
  \item{\dots}{Additional arguments passed to \code{\link{filterVcfBasic}}.}
}

\value{A list with elements \code{vcf}, \code{flag} and \code{flag_comment}. 
\code{vcf} contains the 
filtered \code{CollapsedVCF}, \code{flag} a \code{logical(1)} flag if 
problems were identified, further described in \code{flag_comment}.    }

\author{Markus Riester}



\seealso{\code{\link{filterVcfBasic}}}
\examples{
### This function is typically only called by runAbsolute via the 
### fun.filterVcf and args.filterVcf comments.
library(VariantAnnotation)    
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
vcf <- readVcf(vcf.file, "hg19")
vcf.filtered <- filterVcfMuTect(vcf)        
}
