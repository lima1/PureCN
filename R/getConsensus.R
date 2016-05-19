getConsensus <- function(# Find consensus of multiple PureCN runs
### This function can be used to improve the ranking of PureCN
### purity/ploidy solutions by providing alternate runs. These can
### be either from different samples from the same patient or from
### runAbsoluteCN calls using different parameters.
reference,
### Return object of the runAbsoluteCN() function. 
others
### List of return objects of the runAbsoluteCN() function.
### These are used to order the purity/ploidy solutions of
### the reference object. 
) {
   matching.ids <- lapply(seq_along(reference$results), function(i) {
        ref.obj <- reference$results[[i]]
        sapply(others, function(ret) which.min(sapply(ret$results, 
            function(x) .differenceBySeg(x$seg,ref.obj$seg))))
    })   
    # sum up the likelihoods of the best matches
    sum.ll <- sapply(matching.ids, function(ids) Reduce("+",
        sapply(seq_along(ids), function(j)
        others[[j]]$results[[ids[j]]]$total.log.likelihood)))
    sum.ll <- sum.ll + sapply(reference$results, function(x) 
        x$total.log.likelihood)
    idx <- order(sum.ll, decreasing=TRUE)
    reference.ordered <- reference
    reference.ordered$results <- reference.ordered$results[idx]
    reference.ordered 
}

.differenceBySeg <- function(seg.a, seg.b) {
    gr.a <- GRanges(seqnames=seg.a$chrom, 
            IRanges(start=seg.a$loc.start, end=seg.a$loc.end))
    gr.b <- GRanges(seqnames=seg.b$chrom, 
            IRanges(start=seg.b$loc.start, end=seg.b$loc.end))
    ov <- findOverlaps(gr.a, gr.b)
    mean(abs(seg.a$C[queryHits(ov)] - seg.b$C[subjectHits(ov)]),na.rm=TRUE)
}
