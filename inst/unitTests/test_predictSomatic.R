test_predictSomatic <- function() {
    data(purecn.example.output)
    ret <- predictSomatic(purecn.example.output)
    checkEquals("data.frame", class(ret))
    checkEquals( nrow(
      purecn.example.output$results[[1]]$SNV.posterior$beta.model$posteriors),
                nrow(ret))
}    
