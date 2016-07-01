test_createSNPBlacklist <- function() {
    vcf.file <- system.file("extdata", "example_vcf.vcf", package = "PureCN")
    vcf <- readVcf(vcf.file, "hg19")
    vcf.files <- rep(vcf.file,3)
    ret <- createSNPBlacklist(vcf.files)
    checkEquals("list", class(ret))
    checkEqualsNumeric(start(vcf[rownames(ret[[1]])]), ret[[1]]$start)
    checkEquals(as.character(seqnames(vcf[rownames(ret[[1]])]))
        , as.character(ret[[1]]$chr))
    checkEqualsNumeric(unlist(geno(vcf[rownames(ret[[1]])])$FA[,1]),
        ret[[1]]$Mean.AR)
}    
