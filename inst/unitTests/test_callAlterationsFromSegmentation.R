test_callAlterationsFromSegmentations <- function() {
    data(purecn.example.output)
    seg <- purecn.example.output$results[[1]]$seg
    gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
        package = "PureCN")
    calls <- callAlterationsFromSegmentation(sampleid = seg$ID, 
        chr = seg$chrom, start = seg$loc.start, end = seg$loc.end, 
        num.mark = seg$num.mark, seg.mean = seg$seg.mean, C = seg$C, 
        gc.gene.file = gc.gene.file)

    calls2 <- callAlterations(purecn.example.output)
    checkEquals(sort(rownames(calls2[calls2$type=="AMPLIFICATION",])), 
        sort(rownames(calls$Sample1[calls$Sample1$type=="AMPLIFICATION",])))
}    
