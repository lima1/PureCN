test_getDiploid <- function() {
    data(purecn.example.output)
    ret <- getDiploid(purecn.example.output)
    checkEquals("list", class(ret))
    checkEquals(c("ids", "fraction.non.single"), names(ret))
    checkEquals(1, length(ret$id))
    checkEquals(length(purecn.example.output$results), length(ret$fraction.non.single))
}    
