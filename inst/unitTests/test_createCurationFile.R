test_createCurationFile <- function() {
    data(purecn.example.output)
    file.rds <- 'Sample1_PureCN.rds'
    saveRDS(purecn.example.output, file=file.rds)
    ret <- createCurationFile(file.rds) 
    checkEqualsNumeric( purecn.example.output$results[[1]]$purity , ret$Purity)
    checkEqualsNumeric( purecn.example.output$results[[1]]$ploidy , ret$Ploidy, 
        tolerance=0.1)
    checkTrue(!ret$Curated)
    checkTrue(ret$Flagged)
    checkEquals(purecn.example.output$input$sampleid, as.character(ret$Sampleid))
    
    retx <- readCurationFile(file.rds)
    checkEqualsNumeric( purecn.example.output$results[[1]]$purity , 
        retx$results[[1]]$purity)
    checkEqualsNumeric( purecn.example.output$results[[1]]$ploidy , 
        retx$results[[1]]$ploidy, tolerance=0.1)

    retx <- readCurationFile(file.rds, min.ploidy=2)
    checkEqualsNumeric( purecn.example.output$results[[2]]$purity , 
        retx$results[[1]]$purity)
    checkEqualsNumeric( purecn.example.output$results[[2]]$ploidy , 
        retx$results[[1]]$ploidy, tolerance=0.1)

    retx <- readCurationFile(file.rds, max.ploidy=2)
    checkEquals( rep(TRUE, length(retx$results)) , 
        sapply(retx$results, function(x) x$ploidy) < 2 )

    retx <- readCurationFile(file.rds, report.best.only=TRUE)
    checkEqualsNumeric( purecn.example.output$results[[1]]$purity , 
        retx$results[[1]]$purity)
    checkEqualsNumeric( purecn.example.output$results[[1]]$ploidy , 
        retx$results[[1]]$ploidy, tolerance=0.1)
    checkEqualsNumeric( 1, length(retx$results))
    
    retx <- purecn.example.output
    retx$results[[1]]$purity <- 0.8

    saveRDS(retx, file=file.rds)
    # file already exists
    filename <- file.path(dirname(file.rds), paste(gsub(".rds$", 
        "", basename(file.rds)), "csv", sep = "."))

    createCurationFile(file.rds, overwrite.uncurated=FALSE) 
    ret <- read.csv(filename, as.is=TRUE) 
    checkEqualsNumeric( purecn.example.output$results[[1]]$purity , ret$Purity)
    checkEqualsNumeric( purecn.example.output$results[[1]]$ploidy , ret$Ploidy, 
        tolerance=0.1)

    createCurationFile(file.rds) 
    ret <- read.csv(filename, as.is=TRUE) 
    checkEqualsNumeric( retx$results[[1]]$purity , ret$Purity)
    checkEqualsNumeric( retx$results[[1]]$ploidy , ret$Ploidy, 
        tolerance=0.1)

    ret$Curated <- TRUE
    write.csv(ret, file = filename, row.names = FALSE)
    saveRDS(purecn.example.output, file=file.rds)
    createCurationFile(file.rds)
    ret <- read.csv(filename, as.is=TRUE) 
    checkTrue(ret$Curated)
    checkEqualsNumeric(0.8, ret$Purity)
}    
