test_createCurationFile <- function() {
    data(purecn.example.output)
    ret <- bootstrapSolutions(purecn.example.output) 
    checkEqualsNumeric( purecn.example.output$results[[1]]$purity , ret$results[[1]]$purity)
    checkEqualsNumeric( purecn.example.output$results[[1]]$ploidy , ret$results[[1]]$ploidy, tolerance=0.1)
    checkTrue(length(ret$results) < length(purecn.example.output$results))
}    
