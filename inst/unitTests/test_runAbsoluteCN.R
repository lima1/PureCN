test_runAbsoluteCN <- function() {
    gatk.normal.file <- system.file("extdata", "example_normal.txt", 
        package = "PureCN")
    gatk.normal2.file <- system.file("extdata", "example_normal2.txt", 
        package="PureCN") 
    gatk.normal.files <- c(gatk.normal.file, gatk.normal2.file) 
    gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
        package = "PureCN")
    vcf.file <- system.file("extdata", "example_vcf.vcf", package = "PureCN")
    gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
        package = "PureCN")
    seg.file <- system.file("extdata", "example_seg.txt", 
        package = "PureCN")

    data(purecn.example.output)
    exon.weight.file <- "exon_weights.txt"
    createExonWeightFile(gatk.tumor.file, gatk.normal.files, exon.weight.file)

    # run without a VCF
    ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
        gatk.tumor.file=gatk.tumor.file, 
        candidates=purecn.example.output$candidates, 
        args.segmentation=list(exon.weight.file=exon.weight.file), 
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
    checkException(runAbsoluteCN(gatk.tumor.file, gatk.tumor.file, genome="hg19"),
        silent=TRUE)

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
    checkException(callAlterations(ret))

    # test with numeric chromosome names
    .strip.chr.name <- function(ls) {
        chr.hash <- NULL
        data(chr.hash, envir = environment())
        x <- chr.hash[as.character(ls), 2]
        x[is.na(x)] <- as.numeric(ls[is.na(x)])
        x[is.na(x)] <- max(x, na.rm=TRUE)+1
        x
    }
    vcf <- readVcf(vcf.file, "hg19") 
    vcf <- renameSeqlevels(vcf,as.character( .strip.chr.name(seqlevels(vcf))))
    normCov <- readCoverageGatk( gatk.normal.file )
    tumorCov <- readCoverageGatk( gatk.tumor.file )
    normCov$chr <- .strip.chr.name(normCov$chr)
    tumorCov$chr <- .strip.chr.name(tumorCov$chr)
    normCov$probe <- paste(normCov$chr, ":", normCov$probe_start, "-", normCov$probe_end, sep="")
    tumorCov$probe <- paste(tumorCov$chr, ":", tumorCov$probe_start, "-", tumorCov$probe_end, sep="")

    checkException(runAbsoluteCN(gatk.normal.file=normCov, 
        gatk.tumor.file=tumorCov, remove.off.target.snvs=TRUE,
        gc.gene.file=gc.gene.file,
        vcf.file=vcf, genome="hg19", test.purity=seq(0.3,0.7, by=0.01),
        max.candidate.solutions=2))

    gc3 <- read.delim(gc.gene.file, as.is=TRUE)
    gc3[,1] <- tumorCov$probe
    write.table(gc3, file="tmp3.gc", row.names=FALSE, sep="\t", quote=FALSE)

    ret <- runAbsoluteCN(gatk.normal.file=normCov, 
        gatk.tumor.file=tumorCov, remove.off.target.snvs=TRUE,
        gc.gene.file="tmp3.gc",
        vcf.file=vcf, genome="hg19", test.purity=seq(0.3,0.7, by=0.01),
        max.candidate.solutions=2)

    gpnmb <- callAlterations(ret)["GPNMB",]
    checkEquals("7", gpnmb$chr)
    checkTrue(gpnmb$start>23000000)
    checkTrue(gpnmb$end  <23400000)
    checkTrue(gpnmb$C >= 6) 
    checkEquals("AMPLIFICATION", gpnmb$type) 
    # run the plot function to catch crashes
    plotAbs(ret, 1, type="all")
    checkException(plotAbs(ret, 1, type="BAF", chr="chr1"))
    plotAbs(ret, 1, type="BAF", chr="1")

    # test with a seg.file
    ret <- runAbsoluteCN( gatk.tumor.file = gatk.tumor.file, seg.file=seg.file,
        vcf.file=vcf.file, max.candidate.solutions=1,genome="hg19", 
        test.purity=seq(0.3,0.7, by=0.01))
    
    checkEqualsNumeric(ret$results[[1]]$purity, 0.65, tolerance=0.1)
}    
