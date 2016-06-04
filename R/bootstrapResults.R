.bootstrapResults <- function(results, n=500, top=2, keep=c()) {
    .bootstrapResult <- function(result) {
        lliks <- log(apply(result$SNV.posterior$beta.model$likelihoods[
            !result$SNV.posterior$beta.model$llik.ignored,],1,max))
        lliks <- sum(sample(lliks, replace=TRUE))
        result$log.likelihood + sum(lliks) - 
            sum(result$SNV.posterior$beta.model$llik.ignored)
    }
    best <- lapply(1:n, function(i) head(
        order(sapply(results, .bootstrapResult), decreasing=TRUE), 2))
    bootstrap.value <- sapply(seq_along(results), function(i) 
        sum(sapply(best, function(x) x[1]==i)))/length(best)

    for (i in seq_along(results)) {
        results[[i]]$bootstrap.value <- bootstrap.value[i]
    }    
    best <- unlist(best)
    best <- c(best, keep)
    results <- results[unique(best)]
    results <- results[order(sapply(results, function(x) x$bootstrap.value),decreasing=TRUE)]
#    results[as.numeric(names(sort(table(best),decreasing=TRUE)))]
    results
}

bootstrapResults <- structure(
function(# Filter unlikely purity/ploidy solutions
### This function bootstraps SNVs, then re-ranks solutions 
### by using the bootstrap estimate of the likelihood score, and then keeps 
### only solutions that were ranked highest in any bootstrap replicate.
### Large-scale copy number artifacts can cause true purity/ploidy 
### solutions rank low.
ret,
### Return object of the runAbsoluteCN() function.
n=500,
### Number of bootstrap replicates.
top=2,
### Include solution if it appears in the top n solutions of
### any bootstrap replicate.
keep=c()
### Ids of solutions that should not be removed, even if bootstrap
### value is low, for example diploid solutions. So this will not remove 
### the ret$results[keep] solutions.
 ) {
    if (class(keep) == "logical") keep <- which(keep)
    r <- .bootstrapResults(ret$results, n=n, top=top, keep)
    ret$results <- r
    ret
### Returns the runAbsoluteCN object with low likelihood solutions
### removed. Also adds a bootstrap value to each solution. This value is
### the fraction of bootstrap replicates in which the solution ranked
### first.                
}, ex=function() {
data(purecn.example.output)
ret.boot <- bootstrapResults(purecn.example.output)
plotAbs(ret.boot, type="overview")
})    
