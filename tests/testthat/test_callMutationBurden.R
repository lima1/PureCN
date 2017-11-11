context("callMutationBurden")

data(purecn.example.output)
callableBed <- import(system.file("extdata", "example_callable.bed.gz", 
    package = "PureCN"))

test_that("Example is called correctly", {
    calls <- callMutationBurden(purecn.example.output)
    expect_true(is.na(calls$callable.bases.ontarget))
    exclude <- GRanges(seqnames = "chr1", IRanges(start = 1, 
        end = max(end(callableBed))))
    myVcfFilter <- function(vcf) seqnames(vcf) != "chr2"
    callsCallable <- callMutationBurden(purecn.example.output, 
        callable = callableBed, exclude = exclude, fun.countMutation = myVcfFilter)
    expect_true(callsCallable$callable.bases.ontarget > 0)
    expect_true(callsCallable$callable.bases.flanking > callsCallable$callable.bases.ontarget)
    expect_true(callsCallable$callable.bases.all > callsCallable$callable.bases.flanking)
})    

test_that("Exceptions happen with wrong input", {
    expect_error(callMutationBurden(purecn.example.output, callable = callableBed, 
        exclude = exclude, fun.countMutation = "helloworld"))
    expect_error(callMutationBurden(purecn.example.output, callable = callableBed, 
        exclude = "helloworld"))
    expect_error(callMutationBurden(purecn.example.output, callable = "helloworld"))
})
