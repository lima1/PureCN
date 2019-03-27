#' Call Chromosomal Instability
#' 
#' This function provides detailed CIN information.
#' 
#' 
#' @param res Return object of the \code{\link{runAbsoluteCN}} function.
#' @param id Candidate solution to extract CIN from. \code{id=1} will use the
#' maximum likelihood solution.
#' @param allele.specific Use allele-specific or only total copy number for
#' detecting abnormal regions. Copy-number neutral LOH would be ignored when
#' this parameter is set to \code{FALSE}.
#' @param reference.state Copy number regions different from the reference
#' state are counted as abnormal. Default is \code{dominant} means the most
#' common state. The other option is \code{normal}, which defines normal
#' heterozygous, diploid as reference. The default is robust to errors in
#' ploidy.
#' @return Returns \code{double(1)} with CIN value.
#' @author Markus Riester
#' @seealso \code{\link{runAbsoluteCN}}
#' @examples
#' 
#' data(purecn.example.output)
#' head(callCIN(purecn.example.output))
#' 
#' @export callCIN
callCIN <- function(res, id = 1, allele.specific = TRUE, reference.state = 
                    c("dominant", "normal")) {
    loh <- callLOH(res, id)
    loh$size <- loh$end - loh$start + 1
    # should not happen
    loh <- loh[!is.na(loh$size), ]
    if (allele.specific) loh <- loh[!is.na(loh$M), ]
    reference.state <- match.arg(reference.state)
    loh$state <- if (allele.specific) paste0(loh$C, "/", loh$M) else loh$C
    dominant.state <-  sort(sapply(split(loh$size, loh$state), sum),
                            decreasing=TRUE)[1]
    reference.state.cn <- names(dominant.state)
    if (reference.state == "normal") {
        reference.state.cn <- if (allele.specific) "2/1" else "2"
    }
    loh$is.reference <- loh$state == reference.state.cn
    sum(loh$size[loh$is.reference]) / sum(loh$size)
}
