test_createCurationFile <- function() {
    data(purecn.example.output)
    file.rds <- 'Sample1_PureCN.rds'
    saveRDS(purecn.example.output, file=file.rds)
    ret <- createCurationFile(file.rds) 
    checkEqualsNumeric( purecn.example.output$results[[1]]$purity , ret$Purity)
    checkEqualsNumeric( purecn.example.output$results[[1]]$ploidy , ret$Ploidy, tolerance=0.1)
    checkTrue(!ret$Curated)
    checkTrue(ret$Flagged)
    checkEquals(purecn.example.output$input$sampleid, as.character(ret$Sampleid))
}    
