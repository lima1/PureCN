test_callLOH <- function() {
    data(purecn.example.output)
    ret <- callLOH(purecn.example.output)
    checkEquals("data.frame", class(ret))
    checkEqualsNumeric(7, ncol(ret))
}    
