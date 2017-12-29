context("callAlterationsFromSegmentation")

test_that("Example is called correctly", {
    data(purecn.example.output)
    seg <- purecn.example.output$results[[1]]$seg
    interval.file <- system.file("extdata", "example_intervals.txt", 
        package = "PureCN")
    calls <- callAlterationsFromSegmentation(sampleid = seg$ID, 
        chr = seg$chrom, start = seg$loc.start, end = seg$loc.end, 
        num.mark = seg$num.mark, seg.mean = seg$seg.mean, C = seg$C, 
        interval.file = interval.file)
    calls2 <- callAlterations(purecn.example.output)
    expect_equal(sort(rownames(calls$Sample1[calls$Sample1$type == 
        "AMPLIFICATION", ])), sort(rownames(calls2[calls2$type == 
        "AMPLIFICATION", ])))
})
