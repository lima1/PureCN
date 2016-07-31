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
    cosmic.vcf.file <- system.file("extdata", "example_cosmic.vcf.gz", 
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
    checkException(callLOH(ret))
    checkTrue(grepl("runAbsoluteCN was run without a VCF file", 
        geterrmessage()))
    checkException(runAbsoluteCN(gatk.tumor.file=gatk.tumor.file, 
        genome="hg19"))
    checkTrue(grepl("Need a normal coverage file if log.ratio and seg.file", 
        geterrmessage()))
    checkException(runAbsoluteCN(gatk.tumor.file=gatk.tumor.file, 
        min.ploidy=0, genome="hg19"))
    checkTrue(grepl("min.ploidy or max.ploidy not within expected range", 
        geterrmessage()))

    checkException(runAbsoluteCN(gatk.tumor.file=gatk.tumor.file, 
        max.ploidy=0, genome="hg19"))
    checkTrue(grepl("min.ploidy or max.ploidy not within expected range", 
        geterrmessage()))

    checkException(runAbsoluteCN(gatk.tumor.file=gatk.tumor.file,
        gatk.normal.file=gatk.normal.file, max.ploidy="a", genome="hg19"))
    checkException(runAbsoluteCN(gatk.tumor.file=gatk.tumor.file,
        gatk.normal.file=gatk.normal.file, max.ploidy="a", genome="hg19"))
    checkException(runAbsoluteCN(gatk.tumor.file=gatk.tumor.file,
        gatk.normal.file=gatk.normal.file, test.purity="a", genome="hg19"))
    checkException(runAbsoluteCN(gatk.tumor.file, gatk.tumor.file, 
        genome="hg19"))
    checkTrue(grepl("Tumor and normal are identical", 
        geterrmessage()))
    checkException(runAbsoluteCN(gatk.tumor.file=gatk.tumor.file,
        log.ratio=head(purecn.example.output$input$log.ratio[,2]), 
        genome="hg19"))
    checkTrue(grepl("Length of log.ratio different from tumor coverage",
        geterrmessage()))

    checkException(runAbsoluteCN(gatk.tumor.file=gatk.tumor.file,
        log.ratio=purecn.example.output$input$log.ratio[,2], 
        seg.file=seg.file,
        genome="hg19"))
    checkTrue(grepl("Provide either log.ratio or seg.file, not both",
        geterrmessage()))
    checkException(runAbsoluteCN(gatk.normal.file, gatk.tumor.file,genome="hg19",
        max.homozygous.loss=1.1))
    checkTrue(grepl("max.homozygous.loss not within expected range",
        geterrmessage()))
    checkException(runAbsoluteCN(gatk.normal.file, gatk.tumor.file,genome="hg19",
        prior.K=-1.1))
    checkTrue(grepl("prior.K not within expected range",
        geterrmessage()))
    checkException(runAbsoluteCN(gatk.normal.file, gatk.tumor.file,genome="hg19",
        prior.contamination=c(0.1,-1.1)))
    checkTrue(grepl("prior.contamination not within expected range",
        geterrmessage()))
    
    normalCov <- readCoverageGatk(gatk.normal.file)
    checkException(runAbsoluteCN(normalCov[sample(nrow(normalCov)),], 
        gatk.tumor.file, genome="hg19"))
    checkTrue(grepl("Interval files in normal and tumor different",
        geterrmessage()))
    
    # run with a VCF
    ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
        gatk.tumor.file=gatk.tumor.file, remove.off.target.snvs=TRUE,
        vcf.file=vcf.file, genome="hg19", test.purity=seq(0.3,0.7, by=0.01),
        max.candidate.solutions=2)

    checkEqualsNumeric(ret$results[[1]]$purity, 0.65, tolerance=0.1)
    checkEqualsNumeric(seq(0.3,0.7,by=1/30),
        as.numeric(colnames(ret$candidates$all)))

    vcf <- readVcf(vcf.file, "hg19", param=ScanVcfParam(samples="LIB-02240e4"))

    ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
        gatk.tumor.file=gatk.tumor.file, remove.off.target.snvs=FALSE,
        vcf.file=vcf, genome="hg19", test.purity=seq(0.3,0.7, by=0.01),
        max.candidate.solutions=1)

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
        max.candidate.solutions=1)
    checkTrue(is.na( ret$results[[1]]$gene.calls ))
    checkException(callAlterations(ret))

    # test that correct exons were filtered
    tumor <- readCoverageGatk(gatk.tumor.file)
    log.ratio <- ret$input$log.ratio
    filtered <- cbind(tumor, gc2)[!as.character(tumor$probe) %in% log.ratio$probe,]
    checkTrue(!sum(!(
            filtered$average.coverage < 15 | 
            filtered$targeted.base < 4 | 
            filtered$gc_bias < 0.25 | 
            filtered$gc_bias>0.8)
    ))

    vcf <- readVcf(vcf.file, "hg19") 
    seqlevelsStyle(vcf) <- "ENSEMBL"
    normCov <- readCoverageGatk( gatk.normal.file )
    tumorCov <- readCoverageGatk( gatk.tumor.file )
    normCov$chr <- as.character(normCov$chr) 
    tumorCov$chr <- as.character(tumorCov$chr) 
    seqlevelsStyle(normCov$chr) <- "ENSEMBL"
    seqlevelsStyle(tumorCov$chr) <- "ENSEMBL"
    normCov$probe <- paste(normCov$chr, ":", normCov$probe_start, "-", 
        normCov$probe_end, sep="")
    tumorCov$probe <- paste(tumorCov$chr, ":", tumorCov$probe_start, "-",
        tumorCov$probe_end, sep="")
    
    checkException(runAbsoluteCN(gatk.normal.file, gatk.tumor.file, 
        genome="hg19", vcf.file=vcf))
    checkTrue(grepl("Different chromosome names in coverage and VCF",
        geterrmessage()))

    checkException(runAbsoluteCN(gatk.normal.file=normCov, 
        gatk.tumor.file=tumorCov, remove.off.target.snvs=TRUE,
        gc.gene.file=gc.gene.file,
        vcf.file=vcf, genome="hg19", test.purity=seq(0.3,0.7, by=0.01),
        max.candidate.solutions=2))
    checkTrue(grepl("Intervals of gatk.tumor.file and gc.gene.file do not",
        geterrmessage()))

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

    ret <- runAbsoluteCN( seg.file=seg.file, gc.gene.file=gc.gene.file,
        vcf.file=vcf.file, max.candidate.solutions=1,genome="hg19", 
        test.purity=seq(0.3,0.7, by=0.01),verbose=FALSE)
    
    checkEqualsNumeric(ret$results[[1]]$purity, 0.65, tolerance=0.1)
    tmp <- read.delim(seg.file,as.is=TRUE)
    colnames(tmp)[1:4] <- c("Name", "Chromosome", "Start", "End")
    write.table(tmp, file="seg_wrong.tmp", quote=FALSE, row.names=FALSE,
        sep="\t")
    checkException( runAbsoluteCN( seg.file="seg_wrong.tmp", 
        gc.gene.file=gc.gene.file,
        vcf.file=vcf.file, max.candidate.solutions=1,genome="hg19", 
        test.purity=seq(0.3,0.7, by=0.01),verbose=FALSE))
    
    
    vcf <- readVcf(vcf.file, "hg19")
    # example vcf contains simululated allelic ratios, fix the AD column
    gt1 <- round(geno(vcf)$DP[,1]* as.numeric(geno(vcf)$FA[,1]))
    ad1 <- lapply(seq_along(gt1), function(i) 
        as.numeric(c(geno(vcf)$DP[i,1]-gt1[i], gt1[i])))
    names(ad1) <- names(geno(vcf)$DP[,1])
    geno(vcf)$AD[,1] <- ad1
    geno(vcf)$FA <- NULL
    geno(vcf)$DP <- NULL
    ret <- runAbsoluteCN( gatk.normal.file=gatk.normal.file,
        gatk.tumor.file=gatk.tumor.file, sampleid="LIB-02252e4",
        vcf.file=vcf, max.candidate.solutions=1,genome="hg19", 
        cosmic.vcf.file=cosmic.vcf.file,
        test.purity=seq(0.3,0.7, by=0.01))
    checkEqualsNumeric(ret$results[[1]]$purity, 0.65, tolerance=0.1)
}    
