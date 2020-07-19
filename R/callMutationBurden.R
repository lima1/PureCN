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
#' @param max.prior.somatic Exclude variants with somatic prior
#' probability higher than this cutoff. This is useful for removing
#' hotspot mutations in small panels that might inflate the mutation
#' burden.
#' @param min.cellfraction Exclude variants with cellular fraction
#' lower than this cutoff. These are sub-clonal mutations or artifacts
#' with very low allelic fraction.
#' @param fun.countMutation Function that can be used to filter the
#' input VCF further for filtering, for example to only keep missense
#' mutations. Expects a \code{logical} vector indicating whether variant
#' should be counted (\code{TRUE}) or not (\code{FALSE}). Default
#' is to keep only single nucleotide variants.
#' @param callable \code{GRanges} object with callable genomic regions,
#' for example obtained by \sQuote{GATK CallableLoci} BED file, imported 
#' with \code{rtracklayer}.
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
#' # To calculate exact mutations per megabase, we can provide a BED
#' # file containing all callable regions
#' callableBed <- import(system.file("extdata", "example_callable.bed.gz", 
#'     package = "PureCN"))
#'
#' # We can exclude some regions for mutation burden calculation, 
#' # for example intronic regions. 
#' exclude <- GRanges(seqnames="chr1", IRanges(start=1, 
#'     end=max(end(callableBed))))
#' 
#' # We can also exclude specific mutations by filtering the input VCF
#' myVcfFilter <- function(vcf) seqnames(vcf)!="chr2"
#'
#' callsCallable <- callMutationBurden(purecn.example.output, 
#'     callable=callableBed, exclude=exclude, fun.countMutation=myVcfFilter)
#' 
#' @export callMutationBurden
#' @importFrom stats binom.test
callMutationBurden <- function(res, id = 1, remove.flagged = TRUE,
    min.prior.somatic = 0.1, max.prior.somatic = 1, min.cellfraction = 0,
    fun.countMutation=function(vcf) width(vcf)==1,
    callable = NULL, exclude = NULL) {

    if (is.null(res$input$vcf)) {
        .stopUserError("runAbsoluteCN was run without a VCF file.")
    }

    p <- GRanges(predictSomatic(res, id))

    callableBases <- NA
    callableBasesOntarget <- NA
    callableBasesFlanking <- NA
    if (is.null(callable)) callable <- .estimateCallableRegions(res)
    # calculate the callable genomic region for # mutations/MB calculation
    if (!is.null(callable)) {
        if (!is(callable, "GRanges")) {
            .stopUserError("callable not a GRanges object.")
        } 
        if (!is.null(exclude)) {
            if (!is(exclude, "GRanges")) {
                .stopUserError("exclude not a GRanges object.")
            } 
            callable <- setdiff(callable, exclude)
        }

        targetGRanges <- res$input$log.ratio
        # support for older RDS files
        if (is.null(targetGRanges$on.target)) targetGRanges$on.target <- TRUE
        targetGRanges <- targetGRanges[targetGRanges$on.target]

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
    if (!is.null(fun.countMutation)) {
        if (!is(fun.countMutation, "function")) {
            .stopUserError("fun.countMutation not a function.")        
        }
        vcf <-  res$input$vcf
        vcf <- vcf[which(fun.countMutation(vcf))]
        p <- p[overlapsAny(p, vcf, type = "equal")]
    } 
        
    p <- p[p$prior.somatic >= min.prior.somatic & 
           p$prior.somatic <= max.prior.somatic]
           
    if (remove.flagged) p <- p[!p$FLAGGED]
    
    ret <- data.frame(
        somatic.ontarget=sum(p$ML.SOMATIC & p$on.target==1 & 
            p$CELLFRACTION>min.cellfraction),
        somatic.all=sum(p$ML.SOMATIC & p$CELLFRACTION>min.cellfraction), 
        private.germline.ontarget=sum(!p$ML.SOMATIC & p$on.target==1), 
        private.germline.all=sum(!p$ML.SOMATIC),
        callable.bases.ontarget=callableBasesOntarget,
        callable.bases.flanking=callableBasesFlanking,
        callable.bases.all=callableBases
    )
    
    if (!is.na(ret$callable.bases.ontarget) && ret$callable.bases.ontarget > 0) {
        delta <- ret$callable.bases.ontarget/1e+6
        ci <- binom.test(x = ret$somatic.ontarget, 
            n = ret$callable.bases.ontarget)$conf.int
        ret <- cbind(ret, data.frame(
            somatic.rate.ontarget = ret$somatic.ontarget/delta,
            somatic.rate.ontarget.95.lower = ci[1] * 1e+6,
            somatic.rate.ontarget.95.upper = ci[2] * 1e+6
        ))
        # if multiple CI methods were provided, then ret contains multiple rows
        # now
        ci <- binom.test(x = ret$private.germline.ontarget, 
            n = ret$callable.bases.ontarget)$conf.int
        ret <- cbind(ret, data.frame(
            private.germline.rate.ontarget = ret$private.germline.ontarget/delta,
            private.germline.rate.ontarget.95.lower = ci[1] * 1e+6,
            private.germline.rate.ontarget.95.upper = ci[2] * 1e+6
        ))
    }  
          
    ret
}

.padGranges <- function(target.granges, interval.padding) {
    target.granges.padding <- target.granges
    start(target.granges.padding) <- start(target.granges.padding)-interval.padding
    end(target.granges.padding) <- end(target.granges.padding)+interval.padding
    return(target.granges.padding)
}

.estimateCallableRegions <- function(rds) {
    callable <- NULL
    flog.info("Estimating callable regions.")
    if (is(rds$input$log.ratio, "GRanges")) {
        callable <- reduce(rds$input$log.ratio[rds$input$log.ratio$on.target])
    }    
    return(callable)
}    
