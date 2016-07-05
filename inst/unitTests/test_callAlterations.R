test_callAlterations <- function() {
    data(purecn.example.output)
    calls <- callAlterations(purecn.example.output)
    checkTrue(sum(calls$C < 6 & calls$C > 0.5)==0)
    calls <- callAlterations(purecn.example.output, failed=TRUE)
    checkTrue(sum(calls$gene.mean < 0.9 & calls$gene.mean > -0.9)==0)
    esr2 <- callAlterations(purecn.example.output, all.genes=TRUE)["ESR2",]
    checkEquals("chr14", as.character(esr2$chr))
    checkTrue(esr2$start > 64694600)
    checkTrue(esr2$end < 64761128)
}    
