getDiploid <- structure(function(# Function to extract diploid solutions.
### This function can be used to extract purity and ploidy solutions that
### are diploid with only few CNVs. Since high ploidy solutions 
### typically have a very different copy number profile, one of the 
### identified diploid solutions is likely correct if there are any. This 
### function can be used for automated curation workflows; very silent genomes
### have by definition only a small number of CNVs, making it difficult for the
### algorithm to correctly identify purity and ploidy.
res, 
### Return object of the runAbsoluteCN() function.
min.diploid=0.65, 
### Minimum fraction of genome with normal copy number 2.
max.non.single.gain.loss=0.10
### Maximum fraction of genome with copy number smaller 1 or more than 3.
) {
    cs <- sapply(0:7, function(i) sapply(res$results, function(x) 
                sum(x$seg$size[x$seg$C == i])/sum(x$seg$size)))

    fraction.non.single <- apply(cs[, -(2:4)],1,sum)

    idx <- cs[,3] >= min.diploid & fraction.non.single < max.non.single.gain.loss
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
    diploid <- getDiploid(res)
    if (verbose) message("Found ", length(diploid$ids), " diploid solutions.")

    if (bootstrap && is.null(res$results[[1]]$bootstrap.value)) {
        if (verbose) message("Bootstrapping VCF ",
                "to reduce number of solutions.")
        res <- bootstrapResults(res, keep=diploid$ids)
        diploid <- getDiploid(res)
    }

    ids <- diploid$ids
    # Similar bootstrap value? Then re-order by ploidy closest to diploid
    if (bootstrap && length(diploid$ids)>1) {
        # re-order now, and undo if this looks wrong
        ids <- diploid$ids[order(diploid$fraction.non.single[diploid$ids])]
        if (res$results[[ids[2]]]$bootstrap.value - 
            res$results[[ids[1]]]$bootstrap.value > 0.1 ||
            abs(res$results[[ids[1]]]$purity - 
                res$results[[1]]$purity) > 0.2) {
            ids <- diploid$ids
        }    
    }           
    if (length(ids)>0) {
        res$results <- c(res$results[ids], res$results[-ids])        
    }
    res    
}, ex=function() {
})        
