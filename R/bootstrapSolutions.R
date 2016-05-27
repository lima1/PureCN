.bootstrapSolutions <- function(results, n=500, top=2) {
    .bootstrapSolution <- function(result) {
        lliks <- log(rowMax(result$SNV.posterior$beta.model$likelihoods[!result$SNV.posterior$beta.model$llik.ignored,]))
        lliks <- sum(sample(lliks, replace=TRUE))
        result$log.likelihood + sum(lliks) - sum(result$SNV.posterior$beta.model$llik.ignored)
    }
    best <- lapply(1:n, function(i) head(order(sapply(results, .bootstrapSolution), decreasing=TRUE), 2))
    best <- unlist(best)
    results[as.numeric(names(sort(table(best),decreasing=TRUE)))]
}

bootstrapSolutions <- structure(
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
top=2
### Include solution if it appears in the top n solutions of
### any bootstrap replicate.
 ) {
    r <- .bootstrapSolutions(ret$results, n=n, top=top)
    ret$results <- r
    ret
}, ex=function() {
data(purecn.example.output)
ret.boot <- bootstrapSolutions(purecn.example.output)
plotAbs(ret.boot, type="overview")
})    
