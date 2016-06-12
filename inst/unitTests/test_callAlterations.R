test_callAlterations <- function() {
    data(purecn.example.output)
    calls <- callAlterations(purecn.example.output)
    checkTrue(sum(calls$C < 6 & calls$C > 0.5)==0)
    calls <- callAlterations(purecn.example.output, failed=TRUE)
    checkTrue(sum(calls$gene.mean < 0.9 & calls$gene.mean > -0.9)==0)
}    
