\name{filterVcfMuTect}
\alias{filterVcfMuTect}
\title{Filter VCF MuTect}
\description{Function to remove artifacts and low confidence/quality calls from
a MuTect generated VCF file.}
\usage{filterVcfMuTect(vcf, tumor.id.in.vcf = NULL, stats.file = NULL, 
    ignore = c("clustered_read_position", "fstar_tumor_lod", 
        "nearby_gap_events", "poor_mapping_region_alternate_allele_mapq", 
        "poor_mapping_region_mapq0", "possible_contamination", 
        "strand_artifact"), verbose = TRUE, ...)}
\arguments{
  \item{vcf}{VCF object, read in with the readVcf function from the 
VariantAnnotation package.}
  \item{tumor.id.in.vcf}{The tumor id in the VCF file, optional.}
  \item{stats.file}{MuTect stats file}
  \item{ignore}{
}
  \item{verbose}{Verbose output.}
  \item{\dots}{Additional arguments passed to filterVcfBasic}
}

\value{A list with elements vcf, flag and flag_comment. "vcf" contains the 
filtered CollapsedVCF, "flag" a flag if problems were identified, further 
described in "flag_comment".    }

\author{Markus Riester}




\examples{
### This function is typically only called by runAbsolute via the 
### fun.filterVcf and args.filterVcf comments.
library(VariantAnnotation)    
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
vcf <- readVcf(vcf.file, "hg19")
vcf.filtered <- filterVcfMuTect(vcf)        
}
