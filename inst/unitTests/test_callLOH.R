test_callLOH <- function() {
    data(purecn.example.output)
    ret <- callLOH(purecn.example.output)
    checkEquals("data.frame", class(ret))
    checkEqualsNumeric(7, ncol(ret))
    normal.coverage.file <- system.file('extdata', 'example_normal.txt',
                                            package='PureCN')
    tumor.coverage.file <- system.file('extdata', 'example_tumor.txt',
                                           package='PureCN')
    vcf.file <- system.file('extdata', 'example_vcf.vcf',
                                package='PureCN')
    vcf <- readVcf(vcf.file)
    normal <- readCoverageFile(normal.coverage.file)
    tumor <- readCoverageFile(tumor.coverage.file)
    seqlevelsStyle(vcf) <- "NCBI"
    seqlevelsStyle(normal) <- "NCBI"
    seqlevelsStyle(tumor) <- "NCBI"

    ret <-runAbsoluteCN(normal.coverage.file=normal,
         tumor.coverage.file=tumor, genome='hg19', vcf.file=vcf,
         sampleid='Sample1', 
         min.ploidy=1.4, max.ploidy=2.4, 
         test.purity=seq(0.4,0.7,by=0.05), 
         max.candidate.solutions=1, plot=FALSE)
    loh <- callLOH(ret)
    checkEquals(as.character(1:22), unique(loh$chr))
}    
