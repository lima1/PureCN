test_autoCurateResults <- function() {
    data(purecn.example.output)
    set.seed(123)
    ret <- autoCurateResults(purecn.example.output, bootstrap.n=100)
    checkEquals("list", class(ret))
    checkEquals(names(purecn.example.output), names(ret))
    checkEquals(4, length(ret$results) )
    x <- purecn.example.output
    x$results <- x$results[1]
    ret <- autoCurateResults(x, bootstrap.n=100)
    checkEquals(1, length(ret$results) )
}    
