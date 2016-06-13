test_predictSomatic <- function() {
    data(purecn.example.output)
    ret <- predictSomatic(purecn.example.output)
    checkEquals("data.frame", class(ret))
    checkEquals( nrow(
      purecn.example.output$results[[1]]$SNV.posterior$beta.model$posteriors),
                nrow(ret))
    esr2 <- ret[ret$gene.symbol=="ESR2",]
    checkEquals("chr14", as.character(esr2$chr))
    checkTrue(esr2$start > 64699747)
    checkTrue(esr2$end < 64761128)
}    
