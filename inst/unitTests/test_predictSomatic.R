test_predictSomatic <- function() {
    data(purecn.example.output)
    ret <- predictSomatic(purecn.example.output)
    checkEquals("data.frame", class(ret))
    checkEquals( nrow(
      purecn.example.output$results[[1]]$SNV.posterior$posteriors),
                nrow(ret))
    esr2 <- ret[ret$gene.symbol=="ESR2",]
    checkEquals("chr14", as.character(esr2$chr))
    checkTrue(esr2$start > 64699747)
    checkTrue(esr2$end < 64761128)
    ret <- predictSomatic(purecn.example.output)
    ret.vcf <- predictSomatic(purecn.example.output, return.vcf=TRUE)
    checkEqualsNumeric(ret$start, start(ret.vcf))
    checkEqualsNumeric(ret$end, end(ret.vcf))
    checkEquals(as.character(ret$chr), as.character(seqnames(ret.vcf)))
    checkEqualsNumeric(round(ret$SOMATIC.M1, digits=4), info(ret.vcf)$SM1, tol=0.001)
    checkEqualsNumeric(round(ret$GERMLINE.M1, digits=4), info(ret.vcf)$GM1, tol=0.001)
    checkEqualsNumeric(round(ret$POSTERIOR.SOMATIC, digits=4), info(ret.vcf)$PS, tol=0.001)
    checkEquals(ret$gene.symbol, info(ret.vcf)$GS)
    # test that segments with less than 5 variants are flagged
    flagged <- lapply(split(ret$seg.id, ret$M.SEGMENT.FLAGGED), table)
    checkTrue(min(flagged$`FALSE`) >= 5) 
    checkTrue(max(flagged$`TRUE`) < 5) 
}
