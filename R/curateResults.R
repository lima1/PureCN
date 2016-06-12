getDiploid <- structure(function(# Function to extract diploid solutions.
### This function can be used to extract purity and ploidy solutions that
### are diploid with only few CNVs. Since high ploidy solutions 
### typically have a very different copy number profile, one of the 
### identified diploid solutions is likely correct if there are any. This 
### function can be used for automated curation workflows; very silent genomes
### have by definition only a small number of CNVs, making it difficult for the
### algorithm to correctly identify purity and ploidy. If the maximum
### likelihood solution is diploid, it is always returned; all other solutions
### must pass the more stringent criteria as defined in the function arguments.
res, 
### Return object of the runAbsoluteCN() function.
min.diploid=0.5, 
### Minimum fraction of genome with normal copy number 2.
min.single.gain.loss=0.05,
### Minimum fraction of genome with copy number 1 or 3. This makes sure that 
### low purity samples are not confused with quiet samples.
max.non.single.gain.loss=0.10,
### Maximum fraction of genome with copy number smaller 1 or more than 3.
max.loh=0.5,
### Maximum fraction of genome in LOH. 
min.log.likelihood=NULL
### Minimum copy number log-likelihood to consider sample. If NULL, not tested.
) {
    cs <- sapply(0:7, function(i) sapply(res$results, function(x) 
                sum(x$seg$size[x$seg$C == i])/sum(x$seg$size)))
    
    ploidy <- sapply(res$results, function(x) x$ploidy)
    purity <- sapply(res$results, function(x) x$purity)

    ll <- sapply(res$results, function(x) x$log.likelihood)
    if (is.null(min.log.likelihood)) min.log.likelihood <- min(ll)

    fraction.loh <- sapply(res$results, .getFractionLoh)
    fraction.single <- apply(cs[,  c(2,4)],1,sum)
    fraction.non.single <- apply(cs[, -(2:4)],1,sum)
    
    idx <- ( cs[,3] >= min.diploid & 
        fraction.single >= min.single.gain.loss &
        fraction.non.single < max.non.single.gain.loss & 
        fraction.loh <= max.loh &
        ll >= min.log.likelihood ) |
        ( ploidy > 1.5 & ploidy < 2.5 & seq_along(ploidy)==1 ) 

    ##value<< A list with elements
    list(
        ids=which(idx), ##<< The ids of diploid solutions (res$results[ids]).
        fraction.non.single=fraction.non.single ##<< The fraction of the genome with copy number <1 or >3.
    )    
##end<<    
}, ex=function() {
data(purecn.example.output)
# no diploid solutions in the example
idx <- getDiploid(purecn.example.output)
})    

autoCurateResults <- structure(
function(# Heuristics to find the best purity/ploidy solution.
### This implements a workflow with various heuristics, with the goal
### of identifying correct purity/ploidy solutions in difficult samples.
### This is mainly for automated copy number calling. This function 
### may evolve over time and might produce different rankings
### after PureCN updates.
res,
### Return object of the runAbsoluteCN() function.
bootstrap=TRUE,
### Try to reduce the number of local optima by using the
### bootstrapResults function.
verbose=TRUE
### Verbose output.
) {

    if (bootstrap && is.null(res$results[[1]]$bootstrap.value)) {
        if (verbose) message("Bootstrapping VCF ",
                "to reduce number of solutions.")
        res <- bootstrapResults(res)
    }
    diploid <- getDiploid(res)

    if (verbose) message("Found ", length(diploid$ids), " diploid solutions.")
    ids <- diploid$ids

    # Similar bootstrap value? Then re-order by ploidy closest to diploid
    if (bootstrap && length(diploid$ids)>1) {
        # re-order now
        ids <- diploid$ids[order(diploid$fraction.non.single[diploid$ids])]
        
    # undo re-ordering if bootstrap value of now second solution is 
    # very low, or purity is very different from non-diploid ML solution.
        if ( 
            ( res$results[[ids[1]]]$bootstrap.value < 0.05 &&
              res$results[[ids[2]]]$bootstrap.value > 0.5 ) ||
              abs(res$results[[ids[1]]]$purity - 
                res$results[[1]]$purity) > 0.2 ) {
            ids <- diploid$ids
        }    
    }           
    if (length(ids)>0) {
        res$results <- c(res$results[ids], res$results[-ids])        
    }
    res    
}, ex=function() {
data(purecn.example.output)
# no diploid solutions in the example
example.output.curated <- autoCurateResults(purecn.example.output)
})        
