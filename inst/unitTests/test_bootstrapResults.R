test_bootstrapResults <- function() {
    data(purecn.example.output)
    set.seed(123)
    ret <- bootstrapResults(purecn.example.output, n=100) 
    checkEqualsNumeric( purecn.example.output$results[[1]]$purity, 
        ret$results[[1]]$purity)
    checkEqualsNumeric( purecn.example.output$results[[1]]$ploidy, 
        ret$results[[1]]$ploidy, tolerance=0.1)
    checkTrue(length(ret$results) < length(purecn.example.output$results))
    checkTrue(ret$results[[1]]$bootstrap.value >= 0.5)
    checkTrue(ret$results[[2]]$bootstrap.value < 0.5)
    checkTrue(length(ret$results) >= 2)
    ret <- bootstrapResults(purecn.example.output, n=100, top=3) 
    checkTrue(length(ret$results) >= 3)
}    
