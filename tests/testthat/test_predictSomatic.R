context("predictSomatic")

data(purecn.example.output)
ret <- predictSomatic(purecn.example.output)

test_that("Gene symbol annotation matches", {
    expect_equal(class(ret), "data.frame")
    expect_equal(nrow(ret), nrow(purecn.example.output$results[[1]]$SNV.posterior$posteriors))
    esr2 <- ret[which(ret$gene.symbol == "ESR2"), ]
    expect_equal(as.character(esr2$chr), "chr14")
    expect_true(esr2$start > 64699747)
    expect_true(esr2$end < 64761128)
})

test_that("VCF and data.frame provide equivalent results", {
    ret.vcf <- predictSomatic(purecn.example.output, return.vcf = TRUE)
    expect_equal(start(ret.vcf), ret$start)
    expect_equal(end(ret.vcf), ret$end)
    expect_equal(as.character(seqnames(ret.vcf)), as.character(ret$chr))
    expect_equal(info(ret.vcf)$SM1, round(ret$SOMATIC.M1, digits = 4))
    expect_equal(info(ret.vcf)$GM1, round(ret$GERMLINE.M1, digits = 4))
    expect_equal(info(ret.vcf)$PS, round(ret$POSTERIOR.SOMATIC, 
        digits = 4))
    expect_equal(info(ret.vcf)$GS, ret$gene.symbol)
})

test_that("Segments are flagged", {
    flagged <- lapply(split(ret$seg.id, ret$M.SEGMENT.FLAGGED), 
        table)
    expect_true(min(flagged$`FALSE`) >= 5)
    expect_true(max(flagged$`TRUE`) < 5)
    expect_true(min(ret$M.SEGMENT.POSTERIOR) > 0.5)
    expect_equal(max(ret$M.SEGMENT.POSTERIOR), 1)
})
