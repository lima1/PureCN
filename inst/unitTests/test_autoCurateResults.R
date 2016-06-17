test_autoCurateResults <- function() {
    data(purecn.example.output)
    set.seed(123)
    ret <- autoCurateResults(purecn.example.output, bootstrap.n=100)
    checkEquals("list", class(ret))
    checkEquals(names(purecn.example.output), names(ret))
    checkEquals(2, length(ret$results) )
}    
