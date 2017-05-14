test_callMutationBurden <- function() {
    data(purecn.example.output)
    calls <- callMutationBurden(purecn.example.output)

    checkTrue(is.na(calls$callable.bases.ontarget))

    callableBed <- import(system.file("extdata", "example_callable.bed.gz", 
        package = "PureCN"))

    exclude <- GRanges(seqnames="chr1", IRanges(start=1, 
        end=max(end(callableBed))))
    
    myVcfFilter <- function(vcf) seqnames(vcf)!="chr2"

    callsCallable <- callMutationBurden(purecn.example.output, 
        callable=callableBed, exclude=exclude, fun.countMutation=myVcfFilter)

    checkTrue(callsCallable$callable.bases.ontarget>0)
    checkTrue(callsCallable$callable.bases.flanking > callsCallable$callable.bases.ontarget)
    checkTrue(callsCallable$callable.bases.all > callsCallable$callable.bases.flanking)

    # test input checks
    checkException(callMutationBurden(purecn.example.output,
            callable=callableBed, exclude=exclude, 
            fun.countMutation="helloworld"))
    checkException(callMutationBurden(purecn.example.output,
            callable=callableBed, exclude="helloworld"))
    checkException(callMutationBurden(purecn.example.output,
            callable="helloworld"))
}    
