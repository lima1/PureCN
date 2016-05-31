test_findFocal <- function() {
    data(purecn.example.output)
    ret <- findFocal(purecn.example.output$results[[1]]$seg)
    checkEquals("logical", class(ret))
    checkTrue(nrow(purecn.example.output$results[[1]]$seg)==length(ret))
    checkTrue( min(purecn.example.output$results[[1]]$seg[ret,"C"]) >= 6)
}    
