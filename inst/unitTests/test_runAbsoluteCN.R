test_runAbsoluteCN <- function() {
    normal.coverage.file <- system.file("extdata", "example_normal.txt", 
        package = "PureCN")
    normal2.coverage.file <- system.file("extdata", "example_normal2.txt", 
        package="PureCN") 
    normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file) 
    tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
        package = "PureCN")
    vcf.file <- system.file("extdata", "example_vcf.vcf", package = "PureCN")
    gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
        package = "PureCN")
    seg.file <- system.file("extdata", "example_seg.txt", 
        package = "PureCN")
    cosmic.vcf.file <- system.file("extdata", "example_cosmic.vcf.gz", 
        package = "PureCN")

    data(purecn.example.output)
    target.weight.file <- "exon_weights.txt"
    createTargetWeights(tumor.coverage.file, normal.coverage.files, target.weight.file)

    # run without a VCF
    ret <-runAbsoluteCN(normal.coverage.file=normal.coverage.file, 
        tumor.coverage.file=tumor.coverage.file, 
        candidates=purecn.example.output$candidates, 
        args.segmentation=list(target.weight.file=target.weight.file), 
        max.ploidy=4, max.candidate.solutions=1)

    checkEqualsNumeric(ret$results[[1]]$purity, 0.65, tolerance=0.1)
    
    # chromosomes in segmentation ordered numerically, not alphabetically
    chrom <- ret$results[[1]]$seg$chrom
    checkEqualsNumeric(1:22,chrom[!duplicated(chrom)])

    # test a few exceptions
    checkException(callLOH(ret))
    checkTrue(grepl("runAbsoluteCN was run without a VCF file", 
        geterrmessage()))
    checkException(runAbsoluteCN(tumor.coverage.file=tumor.coverage.file, 
        genome="hg19"))
    checkTrue(grepl("Need a normal coverage file if log.ratio and seg.file", 
        geterrmessage()))
    checkException(runAbsoluteCN(tumor.coverage.file=tumor.coverage.file, 
        min.ploidy=0, genome="hg19"))
    checkTrue(grepl("min.ploidy or max.ploidy not within expected range", 
        geterrmessage()))

    checkException(runAbsoluteCN(tumor.coverage.file=tumor.coverage.file, 
        max.ploidy=0, genome="hg19"))
    checkTrue(grepl("min.ploidy or max.ploidy not within expected range", 
        geterrmessage()))

    checkException(runAbsoluteCN(tumor.coverage.file=tumor.coverage.file,
        normal.coverage.file=normal.coverage.file, max.ploidy="a", genome="hg19"))
    checkException(runAbsoluteCN(tumor.coverage.file=tumor.coverage.file,
        normal.coverage.file=normal.coverage.file, max.ploidy="a", genome="hg19"))
    checkException(runAbsoluteCN(tumor.coverage.file=tumor.coverage.file,
        normal.coverage.file=normal.coverage.file, test.purity="a", genome="hg19"))
    checkException(runAbsoluteCN(tumor.coverage.file, tumor.coverage.file, 
        genome="hg19"))
    checkTrue(grepl("Tumor and normal are identical", 
        geterrmessage()))
    checkException(runAbsoluteCN(tumor.coverage.file=tumor.coverage.file,
        log.ratio=head(purecn.example.output$input$log.ratio[,2]), 
        genome="hg19"))
    checkTrue(grepl("Length of log.ratio different from tumor coverage",
        geterrmessage()))
    checkException(runAbsoluteCN(normal.coverage.file, tumor.coverage.file, 
        prior.purity=1, genome="hg19"))
    checkTrue(grepl("prior.purity must have the same", 
        geterrmessage()))
    checkException(runAbsoluteCN(normal.coverage.file, tumor.coverage.file, 
        prior.purity=c(0.6,1.1), test.purity=c(0.2, 0.6), genome="hg19"))
    checkTrue(grepl("prior.purity not within expected range", 
        geterrmessage()))
    checkException(runAbsoluteCN(normal.coverage.file, tumor.coverage.file, 
        prior.purity=c(0.6,0.9), test.purity=c(0.2, 0.6), genome="hg19"))
    checkTrue(grepl("prior.purity must add to 1. Sum is 1.5", 
        geterrmessage()))


    checkException(runAbsoluteCN(normal.coverage.file, tumor.coverage.file,genome="hg19",
        max.homozygous.loss=1.1))
    checkTrue(grepl("max.homozygous.loss not within expected range",
        geterrmessage()))
    checkException(runAbsoluteCN(normal.coverage.file, tumor.coverage.file,genome="hg19",
        prior.K=-1.1))
    checkTrue(grepl("prior.K not within expected range",
        geterrmessage()))
    checkException(runAbsoluteCN(normal.coverage.file, tumor.coverage.file,genome="hg19",
        prior.contamination=c(0.1,-1.1)))
    checkTrue(grepl("prior.contamination not within expected range",
        geterrmessage()))
    checkException(runAbsoluteCN(normal.coverage.file, tumor.coverage.file,genome="hg19",
        iterations=1))
    checkTrue(grepl("Iterations not in the expected range from",
        geterrmessage()))
    checkException(runAbsoluteCN(normal.coverage.file, tumor.coverage.file,genome="hg19",
        iterations=3000))
    checkTrue(grepl("Iterations not in the expected range from",
        geterrmessage()))
    
    normalCov <- readCoverageGatk(normal.coverage.file)
    checkException(runAbsoluteCN(normalCov[sample(nrow(normalCov)),], 
        tumor.coverage.file, genome="hg19"))
    checkTrue(grepl("Interval files in normal and tumor different",
        geterrmessage()))
    
    # run with a VCF
    ret <-runAbsoluteCN(normal.coverage.file=normal.coverage.file, 
        tumor.coverage.file=tumor.coverage.file, remove.off.target.snvs=TRUE,
        vcf.file=vcf.file, genome="hg19", test.purity=seq(0.3,0.7, by=0.05),
        max.candidate.solutions=1)

    s <- predictSomatic(ret)
    checkEqualsNumeric(s$AR.ADJUSTED, s$AR/s$MAPPING.BIAS)

    checkEqualsNumeric(ret$results[[1]]$purity, 0.65, tolerance=0.1)
    checkEqualsNumeric(seq(0.3,0.7,by=1/30),
        as.numeric(colnames(ret$candidates$all)))

    vcf <- readVcf(vcf.file, "hg19", param=ScanVcfParam(samples="LIB-02240e4"))
    
    myMappingBiasTestFun <- function(vcf, ...) {
         tmp <- rep(1, nrow(vcf))
         idx <- as.logical(seqnames(vcf) == "chr9" & start(vcf)== 35811642)
         tmp[idx] <- 0.9
         tmp
    }    
    ret <-runAbsoluteCN(normal.coverage.file=normal.coverage.file, 
        tumor.coverage.file=tumor.coverage.file, remove.off.target.snvs=FALSE,
        vcf.file=vcf, genome="hg19", test.purity=seq(0.3,0.7, by=0.05),
        max.candidate.solutions=1, fun.setMappingBiasVcf=myMappingBiasTestFun)

    checkEqualsNumeric(ret$results[[1]]$purity, 0.65, tolerance=0.1)
    checkEqualsNumeric(seq(0.3,0.7,by=1/30),
        as.numeric(colnames(ret$candidates$all)))
    
    s <- predictSomatic(ret)
    checkEqualsNumeric(s$AR.ADJUSTED, s$AR/s$MAPPING.BIAS)
    checkEqualsNumeric(0.9, s$MAPPING.BIAS[s$chr=="chr9" & s$start==35811642])

    # test with gc.gene.file without symbols
    gc2 <- read.delim(gc.gene.file, as.is=TRUE)[,-3]
    write.table(gc2, file="tmp.gc", row.names=FALSE, sep="\t", quote=FALSE)
    ret <-runAbsoluteCN(normal.coverage.file=normal.coverage.file, 
        gc.gene.file="tmp.gc",
        tumor.coverage.file=tumor.coverage.file, remove.off.target.snvs=TRUE,
        candidates=purecn.example.output$candidates, 
        vcf.file=vcf.file, genome="hg19", test.purity=seq(0.3,0.7, by=0.05),
        max.candidate.solutions=1)
    checkTrue(is.na( ret$results[[1]]$gene.calls ))
    checkException(callAlterations(ret))

    # test that correct exons were filtered
    tumor <- readCoverageGatk(tumor.coverage.file)
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
    normCov <- readCoverageGatk( normal.coverage.file )
    tumorCov <- readCoverageGatk( tumor.coverage.file )
    normCov$chr <- as.character(normCov$chr) 
    tumorCov$chr <- as.character(tumorCov$chr) 
    seqlevelsStyle(normCov$chr) <- "ENSEMBL"
    seqlevelsStyle(tumorCov$chr) <- "ENSEMBL"
    normCov$probe <- paste(normCov$chr, ":", normCov$probe_start, "-", 
        normCov$probe_end, sep="")
    tumorCov$probe <- paste(tumorCov$chr, ":", tumorCov$probe_start, "-",
        tumorCov$probe_end, sep="")
    
    checkException(runAbsoluteCN(normal.coverage.file, tumor.coverage.file, 
        genome="hg19", vcf.file=vcf))
    checkTrue(grepl("Different chromosome names in coverage and VCF",
        geterrmessage()))

    checkException(runAbsoluteCN(normal.coverage.file=normCov, 
        tumor.coverage.file=tumorCov, remove.off.target.snvs=TRUE,
        gc.gene.file=gc.gene.file,
        vcf.file=vcf, genome="hg19", test.purity=seq(0.3,0.7, by=0.05),
        max.candidate.solutions=1))
    checkTrue(grepl("Intervals of tumor.coverage.file and gc.gene.file do not",
        geterrmessage()))

    gc3 <- read.delim(gc.gene.file, as.is=TRUE)
    gc3[,1] <- tumorCov$probe
    write.table(gc3, file="tmp3.gc", row.names=FALSE, sep="\t", quote=FALSE)

    ret <- runAbsoluteCN(normal.coverage.file=normCov, 
        tumor.coverage.file=tumorCov, remove.off.target.snvs=TRUE,
        gc.gene.file="tmp3.gc",
        vcf.file=vcf, genome="hg19", test.purity=seq(0.3,0.7, by=0.05),
        max.candidate.solutions=1)

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
    ret <- runAbsoluteCN( tumor.coverage.file = tumor.coverage.file, seg.file=seg.file,
        vcf.file=vcf.file, max.candidate.solutions=1,genome="hg19", 
        test.purity=seq(0.3,0.7, by=0.05))
    checkEqualsNumeric(ret$results[[1]]$purity, 0.65, tolerance=0.1)

    ret <- runAbsoluteCN( seg.file=seg.file, gc.gene.file=gc.gene.file,
        vcf.file=vcf.file, max.candidate.solutions=1,genome="hg19", 
        test.purity=seq(0.3,0.7, by=0.05),verbose=FALSE)
    
    checkEqualsNumeric(ret$results[[1]]$purity, 0.65, tolerance=0.1)
    tmp <- read.delim(seg.file,as.is=TRUE)
    colnames(tmp)[1:4] <- c("Name", "Chromosome", "Start", "End")
    write.table(tmp, file="seg_wrong.tmp", quote=FALSE, row.names=FALSE,
        sep="\t")
    checkException( runAbsoluteCN( seg.file="seg_wrong.tmp", 
        gc.gene.file=gc.gene.file,
        vcf.file=vcf.file, max.candidate.solutions=1,genome="hg19", 
        test.purity=seq(0.3,0.7, by=0.05),verbose=FALSE))
    
    # test with a log.ratio and no tumor file
    log.ratio <- calculateLogRatio(readCoverageGatk(normal.coverage.file),
        readCoverageGatk(tumor.coverage.file), verbose=FALSE)

    ret <- runAbsoluteCN( log.ratio=log.ratio,
        gc.gene.file=gc.gene.file,
        vcf.file=vcf.file, max.candidate.solutions=1,genome="hg19", 
        test.purity=seq(0.3,0.7, by=0.05))
    checkEqualsNumeric(ret$results[[1]]$purity, 0.65, tolerance=0.1)
    ret <- runAbsoluteCN( log.ratio=log.ratio, seg.file=seg.file,
        gc.gene.file=gc.gene.file,
        vcf.file=vcf.file, max.candidate.solutions=1,genome="hg19", 
        test.purity=seq(0.3,0.7, by=0.05))
    checkEqualsNumeric(ret$results[[1]]$purity, 0.65, tolerance=0.1)
     
    vcf <- readVcf(vcf.file, "hg19")
    # example vcf contains simululated allelic ratios, fix the AD column
    gt1 <- round(geno(vcf)$DP[,1]* as.numeric(geno(vcf)$FA[,1]))
    ad1 <- lapply(seq_along(gt1), function(i) 
        as.numeric(c(geno(vcf)$DP[i,1]-gt1[i], gt1[i])))
    names(ad1) <- names(geno(vcf)$DP[,1])
    geno(vcf)$AD[,1] <- ad1
    geno(vcf)$FA <- NULL
    geno(vcf)$DP <- NULL
    ret <- runAbsoluteCN( normal.coverage.file=normal.coverage.file,
        tumor.coverage.file=tumor.coverage.file, sampleid="LIB-02252e4",
        vcf.file=vcf, max.candidate.solutions=1,genome="hg19", 
        cosmic.vcf.file=cosmic.vcf.file,
        test.purity=seq(0.3,0.7, by=0.01))
    checkEqualsNumeric(ret$results[[1]]$purity, 0.65, tolerance=0.1)
    # test min.ploidy bug
    ret <- runAbsoluteCN(normal.coverage.file, tumor.coverage.file, min.ploidy=2.2, 
        max.ploidy=4, genome="hg19", test.purity=seq(0.3,0.7, by=0.05), 
        plot.cnv=FALSE, max.candidate.solutions=1)

    checkTrue(ret$results[[1]]$ploidy > 2)
    checkTrue(ret$results[[1]]$ploidy < 4)

    # check filterTargets
    checkException(runAbsoluteCN(normal.coverage.file = normal.coverage.file,  
        tumor.coverage.file = tumor.coverage.file, genome = "hg19", 
        args.filterTargets=list(normalDB=vcf.file)))
    checkTrue(grepl("normalDB not a valid normalDB object", 
        geterrmessage()))
    
    normalDB <- createNormalDatabase(normal.coverage.files)

    tmp <- normalDB
    tmp$normal.coverage.files <- NULL
    checkException(runAbsoluteCN(normal.coverage.file = normal.coverage.file,  
        tumor.coverage.file = tumor.coverage.file, genome = "hg19", 
        args.filterTargets=list(normalDB=tmp)))
    checkTrue(grepl("normalDB appears to be empty", 
        geterrmessage()))
    ret <- runAbsoluteCN(normal.coverage.file = normal.coverage.file, 
        tumor.coverage.file = tumor.coverage.file, genome = "hg19", vcf.file = vcf.file, 
        sampleid = "Sample1", gc.gene.file = gc.gene.file, 
        args.filterTargets = list(normalDB = normalDB,
        filter.lowhigh.gc=0
        ), 
        plot.cnv=FALSE,
        max.ploidy = 3, test.purity = seq(0.4, 0.7, by = 0.05), 
        max.candidate.solutions = 1)
    tumor <- readCoverageGatk(tumor.coverage.file)
    idx <- tumor$probe %in% ret$input$log.ratio$probe
    cutoff <- median(normalDB$exon.median.coverage)*0.3

    checkEqualsNumeric(0, sum(!(normalDB$exon.median.coverage[!idx] < cutoff | 
        tumor$targeted.base[!idx] < 4 | tumor$average.coverage[!idx] < 14)))

   checkTrue( sum(!(normalDB$exon.median.coverage[idx] < cutoff | 
    tumor$targeted.base[idx] < 4 | tumor$average.coverage[idx] < 14)) > 9000)
    

    # run with minimal segmentation function:
    testSeg <- function(seg, ...) return(seg)

    res <- runAbsoluteCN(normal.coverage.file, tumor.coverage.file, 
        seg.file=seg.file, fun.segmentation=testSeg, max.ploidy = 4, 
        test.purity = seq(0.3, 0.7, by = 0.05), max.candidate.solutions=1,
        genome='hg19')

    seg <- read.delim(seg.file)
    checkEqualsNumeric(nrow(seg), nrow(res$results[[1]]$seg))
    checkEquals(seg$seg.mean, res$results[[1]]$seg$seg.mean)
}    
