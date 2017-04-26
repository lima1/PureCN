#' Call mutation burden
#' 
#' This function provides detailed mutation burden information.
#' 
#' 
#' @param res Return object of the \code{\link{runAbsoluteCN}} function.
#' @param id Candidate solution to extract mutation burden from. 
#' \code{id=1} will use the maximum likelihood solution.
#' @param remove.flagged Remove variants flagged by 
#' \code{\link{predictSomatic}}.
#' @param min.prior.somatic Exclude variants with somatic prior
#' probability lower than this cutoff.
#' @param min.cellfraction Exclude variants with cellular fraction
#' lower than this cutoff. These are sub-clonal mutations or artifacts
#' with very low allelic fraction.
#' @param fun.countMutation Function that can be used to filter the
#' input VCF further for filtering, for example to only keep missense
#' mutations. Expects a \code{logical} vector indicating whether variant
#' should be counted (\code{TRUE}) or not (\code{FALSE}).
#' @param callable \code{GRanges} object with callable genomic regions,
#' for example obtained by GATK CallableLoci BED file, imported with
#' \code{rtracklayer}.
#' @param exclude \code{GRanges} object with genomic regions that
#' should be excluded from the \code{callable} regions, for example
#' intronic regions. Requires \code{callable}.
#' @return Returns \code{data.frame} with mutation counts and sizes
#' of callable regions.
#' @author Markus Riester
#' @seealso \code{\link{runAbsoluteCN}} \code{\link{predictSomatic}}
#' @examples
#' 
#' data(purecn.example.output)
#' callMutationBurden(purecn.example.output)
#' 
#' @export callMutationBurden
callMutationBurden <- function(res, id = 1, remove.flagged = TRUE, 
    min.prior.somatic=0.1, min.cellfraction=0, fun.countMutation=NULL,
    callable=NULL, exclude=NULL) {

    if (is.null(res$input$vcf)) {
        .stopUserError("runAbsoluteCN was run without a VCF file.")
    }

    p <- GRanges(predictSomatic(res, id))

    callableBases <- NA
    callableBasesOntarget <- NA
    callableBasesFlanking <- NA
    
    # calculate the callable genomic region for # mutations/MB calculation
    if (!is.null(callable)) {
        if (class(callable) != "GRanges") {
            .stopUserError("callable not a GRanges object.")
        } 
        if (!is.null(exclude)) {
            if (class(exclude) != "GRanges") {
                .stopUserError("exclude not a GRanges object.")
            } 
            callable <- setdiff(callable, exclude)
        }

        targetGRanges <- GRanges(res$input$log.ratio$probe)

        intervalPadding <- res$input$args$filterVcf$interval.padding
        if (is.null(intervalPadding)) intervalPadding <- 50
        targetGRangesPadding <- .padGranges(targetGRanges, intervalPadding)

        callableOntarget <- intersect(targetGRanges, callable)
        callableFlanking <- intersect(targetGRangesPadding, callable)
        callableBases <- sum(as.numeric(width(reduce(callable))))
        callableBasesOntarget <- sum(as.numeric(width(reduce(callableOntarget))))
        callableBasesFlanking <- sum(as.numeric(width(reduce(callableFlanking))))
        p <- p[overlapsAny(p, callable)]
    }    
    
    # filter mutations, for example if the user wants to 
    # calculate missense burden
    p$countMutation <- TRUE
    if (!is.null(fun.countMutation)) {
        if (class(fun.countMutation) != "function") {
            .stopUserError("fun.countMutation not a function.")        
        }
        vcf <-  res$input$vcf
        vcf <- vcf[overlapsAny(vcf, p)]
        p$countMutation <- fun.countMutation(vcf)
    } 
        
    p <- p[p$prior.somatic >= min.prior.somatic]
    if (remove.flagged) p <- p[!p$FLAGGED]
        
    data.frame(
        somatic.ontarget=sum(p$ML.SOMATIC & p$on.target==1 & 
            p$CELLFRACTION>min.cellfraction & p$countMutation),
        somatic.all=sum(p$ML.SOMATIC & 
            p$CELLFRACTION>min.cellfraction & p$countMutation), 
        private.germline.ontarget=sum(!p$ML.SOMATIC & p$on.target==1), 
        private.germline.all=sum(!p$ML.SOMATIC),
        callable.bases.ontarget=callableBasesOntarget,
        callable.bases.flanking=callableBasesFlanking,
        callable.bases.all=callableBases
    )    
}

.padGranges <- function(target.granges, interval.padding) {
    target.granges.padding <- target.granges
    start(target.granges.padding) <- start(target.granges.padding)-interval.padding
    end(target.granges.padding) <- end(target.granges.padding)+interval.padding
    return(target.granges.padding)
}

