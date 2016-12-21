#' Find focal amplifications
#' 
#' Function to find focal amplifications in segmented data.  This is
#' automatically called in \code{\link{runAbsoluteCN}}.
#' 
#' 
#' @param seg Segmentation data.
#' @param max.size Cutoff for focal in base pairs.
#' @param cn.diff Minimum copy number delta between neighboring segments.
#' @param min.amp.cn Minimum amplification integer copy number. Segments with
#' lower copy number are not tested.
#' @return \code{logical(n)}, indicating for all n segments whether they are
#' focally amplified or not.
#' @author Markus Riester
#' @seealso \code{\link{runAbsoluteCN}}
#' @examples
#' 
#' normal.coverage.file <- system.file("extdata", "example_normal.txt", 
#'     package="PureCN")
#' tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
#'     package="PureCN")
#' vcf.file <- system.file("extdata", "example_vcf.vcf", 
#'     package="PureCN")
#' gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
#'     package="PureCN")
#' 
#' # The max.candidate.solutions, max.ploidy and test.purity parameters are set to
#' # non-default values to speed-up this example.  This is not a good idea for real
#' # samples.
#' ret <-runAbsoluteCN(normal.coverage.file=normal.coverage.file, 
#'     tumor.coverage.file=tumor.coverage.file, vcf.file=vcf.file, genome="hg19", 
#'     sampleid="Sample1", gc.gene.file=gc.gene.file,
#'     max.candidate.solutions=1, max.ploidy=4, test.purity=seq(0.3,0.7,by=0.05), 
#'     args.focal=list(max.size = 2e+06), fun.focal=findFocal)
#' 
#' @export findFocal
findFocal <- function(seg, max.size = 3000000, cn.diff = 2, min.amp.cn = 5) {
    focal <- rep(FALSE, nrow(seg))
    for (i in seq_len(nrow(seg))) {
        if (seg$C[i] < min.amp.cn) next
        if (seg$size[i] > max.size) next
        size <- seg$size[i]
        if (i>1) {
            for (j in (i-1):1) {
                if (seg$C[j] < seg$C[i]-cn.diff) {
                    break
                }
                size <- size + seg$size[j]    
            }        
        }    
        if (i<nrow(seg)) {
            for (j in (i+1):nrow(seg)) {
                if (seg$C[j] < seg$C[i]-cn.diff) {
                    break
                }
                size <- size + seg$size[j]    
            }        
        }    
        focal[i] <- size < max.size        
    }
    focal
}    
