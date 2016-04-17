\name{createSNPBlacklist}
\alias{createSNPBlacklist}
\title{Create SNP black list}
\description{Function to create a black list of germline SNPs with expected allelic
fraction (AF) smaller than 0.5 in diploid genomes.  }
\usage{createSNPBlacklist(vcf.files, n = min(10, length(vcf.files)), 
    low.af = 0.025, high.af = 0.1, genome = "hg19")}
\arguments{
  \item{vcf.files}{List of VCF files. When a VCF file contains multiple samples, 
it will ignore all samples except the first.}
  \item{n}{Required number of VCF files showing low allelic fraction to 
blacklist a SNP id. }
  \item{low.af}{Defines a low AF p-value.}
  \item{high.af}{Defines a high AF p-value. For every sample with high AF p-value, 
there must be one more sample with low AF to reach the cutoff.  }
  \item{genome}{Version of the reference genome, required for the readVcf() function.}
}

\value{A list with elements snp.black.list and segmented. 
"snp.black.list" is just a list of SNP ids.
"segmented" blacklists whole regions.}

\author{Markus Riester}




\examples{
# Assume VCF files of normals (for example obtained by a MuTect artifact
# detection run) are in directory poolofnormals:
mutect.normal.files <- dir("poolofnormals", pattern="vcf$", full.names=TRUE) 

# These files do not exist in our example, so we do not run the function here.
#snp.blacklist <- createSNPBlacklist(mutect.normal.files)
}
