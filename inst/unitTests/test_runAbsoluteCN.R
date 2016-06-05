test_runAbsoluteCN <- function() {
    gatk.normal.file <- system.file("extdata", "example_normal.txt", 
        package = "PureCN")
    gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
        package = "PureCN")
    vcf.file <- system.file("extdata", "example_vcf.vcf", package = "PureCN")
    gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
        package = "PureCN")
    data(purecn.example.output)

    # run without a VCF
    ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
        gatk.tumor.file=gatk.tumor.file, 
        candidates=purecn.example.output$candidates, 
        max.candidate.solutions=2)

    checkEqualsNumeric(ret$results[[1]]$purity, 0.65, tolerance=0.1)
    
    # chromosomes in segmentation ordered numerically, not alphabetically
    chrom <- ret$results[[1]]$seg$chrom
    checkEqualsNumeric(1:22,chrom[!duplicated(chrom)])

    # test a few exceptions
    checkException(runAbsoluteCN(gatk.tumor.file=gatk.tumor.file), silent=TRUE)
    checkException(runAbsoluteCN(gatk.tumor.file=gatk.tumor.file, min.ploidy=0),
        silent=TRUE)
    checkException(runAbsoluteCN(gatk.tumor.file=gatk.tumor.file, max.ploidy=0),
        silent=TRUE)
    checkException(runAbsoluteCN(gatk.tumor.file=gatk.tumor.file,
        gatk.normal.file=gatk.normal.file, max.ploidy="a"), silent=FALSE)
    checkException(runAbsoluteCN(gatk.tumor.file=gatk.tumor.file,
        gatk.normal.file=gatk.normal.file, max.ploidy="a"), silent=FALSE)
    checkException(runAbsoluteCN(gatk.tumor.file=gatk.tumor.file,
        gatk.normal.file=gatk.normal.file, test.purity="a"), silent=FALSE)

    # run with a VCF
    ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
        gatk.tumor.file=gatk.tumor.file, remove.off.target.snvs=TRUE,
        vcf.file=vcf.file, genome="hg19", test.purity=seq(0.3,0.7, by=0.01),
        max.candidate.solutions=2)

    checkEqualsNumeric(ret$results[[1]]$purity, 0.65, tolerance=0.1)
    checkEqualsNumeric(seq(0.3,0.7,by=1/30),
        as.numeric(colnames(ret$candidates$all)))

    ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
        gatk.tumor.file=gatk.tumor.file, remove.off.target.snvs=FALSE,
        vcf.file=vcf.file, genome="hg19", test.purity=seq(0.3,0.7, by=0.01),
        max.candidate.solutions=2)

    checkEqualsNumeric(ret$results[[1]]$purity, 0.65, tolerance=0.1)
    checkEqualsNumeric(seq(0.3,0.7,by=1/30),
        as.numeric(colnames(ret$candidates$all)))

    # test with gc.gene.file without symbols
    gc2 <- read.delim(gc.gene.file, as.is=TRUE)[,-3]
    write.table(gc2, file="tmp.gc", row.names=FALSE, sep="\t", quote=FALSE)
    ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
        gc.gene.file="tmp.gc",
        gatk.tumor.file=gatk.tumor.file, remove.off.target.snvs=TRUE,
        candidates=purecn.example.output$candidates, 
        vcf.file=vcf.file, genome="hg19", test.purity=seq(0.3,0.7, by=0.01),
        max.candidate.solutions=2)
    checkTrue(is.na( ret$results[[1]]$gene.calls ))
}    
