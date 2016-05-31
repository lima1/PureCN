test_autoCurateResults <- function() {
    data(purecn.example.output)
    ret <- autoCurateResults(purecn.example.output)
    checkEquals("list", class(ret))
    checkEquals(names(purecn.example.output), names(ret))
    checkEquals(2, length(ret$results) )
}    
