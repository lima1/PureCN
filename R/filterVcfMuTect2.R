#' Filter VCF MuTect2
#' 
#' Function to remove artifacts and low confidence/quality calls from a 
#' GATK4/MuTect2 generated VCF file. Also applies filters defined in 
#' \code{filterVcfBasic}.
#' 
#' 
#' @param vcf \code{CollapsedVCF} object, read in with the \code{readVcf}
#' function from the VariantAnnotation package.
#' @param tumor.id.in.vcf The tumor id in the VCF file, optional.
#' @param ignore MuTect2 flags that mark variants for exclusion.
#' @param \dots Additional arguments passed to \code{\link{filterVcfBasic}}.
#' @return A list with elements \code{vcf}, \code{flag} and
#' \code{flag_comment}.  \code{vcf} contains the filtered \code{CollapsedVCF},
#' \code{flag} a \code{logical(1)} flag if problems were identified, further
#' described in \code{flag_comment}.
#' @author Markus Riester
#' @seealso \code{\link{filterVcfBasic}}
#' @examples
#' 
#' ### This function is typically only called by runAbsolute via the 
#' ### fun.filterVcf and args.filterVcf comments.
#' library(VariantAnnotation)    
#' vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
#' vcf <- readVcf(vcf.file, "hg19")
#' vcf.filtered <- filterVcfMuTect(vcf)        
#' 
#' @export filterVcfMuTect2
filterVcfMuTect2 <- function(vcf, tumor.id.in.vcf = NULL,
ignore=c("clustered_events", "t_lod", "str_contraction", 
"read_position", "fragment_length", "multiallelic", "clipping",
"strand_artifact"),
... ){
    if (is.null(fixed(vcf)$FILTER)) return(
        filterVcfBasic(vcf, tumor.id.in.vcf, ...))
    
    n <- nrow(vcf)

    ids <- sort(unique(unlist(sapply(ignore, grep, fixed(vcf)$FILTER))))
    if (length(ids)) vcf <- vcf[-ids]
    flog.info("Removing %i MuTect2 calls due to blacklisted failure reasons.", 
        n-nrow(vcf))
    filterVcfBasic(vcf, tumor.id.in.vcf, ...)
}
