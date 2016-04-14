# Quick Start

## Basic Input Files 

You will need coverage files in GATK DepthOfCoverage format for tumor and a control sample:

    Target  total_coverage  average_coverage    Sample1_total_cvg   Sample1_mean_cvg    Sample1_granular_Q1 Sample1_granular_median Sample1_granular_Q3 Sample1_._above_15
    chr1:69091-70009    0   0   0   0   1   1   1   0
    chr1:367659-368598  6358    9.25121153819759    6358    9.25121153819759    1   7   13  11.6
    chr1:621096-622035  6294    9.16910019318401    6294    9.16910019318401    1   7   12  9.5

PureCN will look for the columns "Target", "total_coverage", and
"average_coverage". The algorithm works best when the coverage data is already
GC normalized. The steps below will NOT GC normalize.

You will also need a VCF file with both somatic and germline variants, for
example as obtained by MuTect.

## First Run

In R:

    library(PureCN)
    ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
                        gatk.tumor.file=gatk.tumor.file, 
                        vcf.file=vcf.file)

To obtain gene level copy number calls, we need an exon-to-gene mapping file.
This is provided together with GC content via the gc.gene.file argument. The
expected format:

    Target	gc_bias	Gene
    chr1:69091-70009	0.427638737758433	OR4F5
    chr1:367659-368598	0.459574468085106	OR4F29
    chr1:621096-622035	0.459574468085106	OR4F3

## Pool of Normals

If a pool of normal samples is available, we can use the coverage data to
estimate the expected variance in coverage per exon. This will help the
segmentation algorithm to filter noise:

    gatk.normal.files <- dir("poolofnormals", 
        pattern="coverage.sample_interval_summary", full.names=TRUE) 
    createExonWeightFile(gatk.normal.files[1:2], gatk.normal.files[-(1:2)], "exon_weights.txt")

This will calculate log-ratios for all normals using the first two as tumor and
the remaining ones as normal. This estimates then exon level coverage standard
deviations.

We can also use the pool of normal to find SNPs with biased allelic fractions
(significantly different from 0.5 for heterozygous SNPs):

    mutect.normal.files <- dir("poolofnormals", pattern="vcf$", full.names=TRUE) 
    snp.blacklist <- createSNPBlacklist(mutect.normal.files)
    write.csv(snp.bl, file="SNP_blacklist.csv")
    write.csv(snp.bl, file="SNP_blacklist_segmented.csv", row.names=FALSE, quote=FALSE)

Finally, we can use the pool of normals to find normal samples for log-ratio
calculation, especially when no matched normal sample is available:

    if (!file.exists("normalDB.rds")) {
        normalDB <- createNormalDatabase(gatk.normal.files)
        saveRDS(normalDB, "normalDB.rds")
    } else {
        normalDB <- readRDS("normalDBe4.rds")
    }
    
    # get the best 4 normals and average them
    gatk.normal.file <- findBestNormal(gatk.tumor.file, normalDB, num.normals=4)
    pool <- poolCoverage(lapply(gatk.normal.file, read.coverage.gatk), remove.chrs=c('chrX', 'chrY'))
    
    # get the best normal
    gatk.normal.file <- findBestNormal(gatk.tumor.file, normalDB)

## VCF Artifact Filtering

If MuTect was run in matched normal mode, then all germline variants are
rejected, that means we cannot just filter by the PASS/REJECT MuTect flags. The
filterVcfMuTect function optionally reads the MuTect stats file and will keep
germline variants, while removing potential artifacts. Otherwise it will use
only the filters based on read depths as defined in filterVcfBasic.

## Recommended Run

Finally, we can run PureCN with all that information:

      ret <-runAbsoluteCN(gatk.normal.file=pool, gatk.tumor.file=gatk.tumor.file, 
        vcf.file=vcf.file, sampleid='Sample1', gc.gene.file=gc.gene.file, 
        args.filterVcf=list(snp.blacklist=snp.blacklist, stats.file=mutect.stats.file), 
        args.segmentation=list(exon.weight.file=exon.weight.file), 
        post.optimize=TRUE)

The post.optimize flag will increase the runtime by a lot, but might be worth
it if the copy number profile is not very clean. This is recommended with lower
coverage datasets or if there are significant capture biases. For high quality
whole exome data, this is typically not necessary.

We also need to create a few output files:

    file.rds <- 'Sample1_PureCN.rds'
    saveRDS(ret, file=file.rds)
    pdf('Sample1_PureCN.pdf', width=10, height=12)
    plotAbs(ret, type='all')
    dev.off()

## Manual Curation

For prediction of SNV status (germline vs somatic, sub-clonal vs. clonal,
homozyous vs heterozygous), it is important that both purity and ploidy are
correct. We provide functionality for curating results:

    createCurationFile(file.rds) 

This will generate a CSV file, in which the correct purity and ploidy values
can be manually entered. It also contains a column "Curated", which should be
set to TRUE, otherwise the file will be overwritten when re-run.

Then in R, the correct solution (closest to the combination in the CSV file)
can be loaded with the readCurationFile function:

    ret <- readCurationFile(file.rds)
    
## Custom Segmentation

By default, we will use DNAcopy to segment the log-ratio. You can easily change
that to your favorite method and the segmentationCBS function can serve as an
example.

    segmentationCBS <- function (normal, tumor, log.ratio, plot.cnv, coverage.cutoff, 
        sampleid = sampleid, exon.weight.file = NULL, alpha = 0.005, 
        vcf = NULL, tumor.id.in.vcf = 1, verbose = TRUE, ...) 

It is also possible to provide segmented file, which we however only recommend
when matched SNP6 data is available (otherwise it is better to customize the
segmentation function as described above). The expected file format is:

    ID  chrom   loc.start   loc.end num.mark    seg.mean
    Sample1   1   61723   5773942 2681    0.125406444072723
    Sample1   1   5774674 5785170 10  -0.756511807441712

## Output

The plotAbs() call above will generate the main plots shown in the manuscript. The R data file (file.rds) contains gene level copy number calls, SNV status and LOH calls.
The purity/ploidy combinations are sorted by likelihood at stored in ret$results.

    names(ret)
    head(ret$results[[1]]$gene.calls, 3)
    chr    start      end C seg.mean seg.id number.exons gene.mean   gene.min
    TP73 chr1  3598833  3649721 3   0.3434      1           14 0.1755614 -0.3861084
    CHD5 chr1  6161904  6240126 3   0.3434      1           44 0.3059672 -0.8812580
    MTOR chr1 11166594 11319542 3   0.3434      1           62 0.2360036 -0.1346406
    gene.max focal num.snps.loh.segment percentage.loh.in.loh.segment
    TP73 0.6437573 FALSE                    0                             0
    CHD5 1.5984124 FALSE                   49                             0
    MTOR 0.6705342 FALSE                   49                             0

This data.frame also contains gene level LOH information. The SNV posteriors:

    head(ret$results[[1]]$SNV.posterior$beta.model$posteriors, 3)
               seqnames   start     end SOMATIC.M0   SOMATIC.M1   SOMATIC.M2
    rs6696489      chr1 6162054 6162054          0 9.155037e-04 2.938312e-07
    rs59788818     chr1 6171992 6171992          0 2.364795e-05 3.397817e-07
    rs2843493      chr1 6184092 6184092          0 9.362524e-16 2.714454e-04
                 SOMATIC.M3   SOMATIC.M4    SOMATIC.M5    SOMATIC.M6 SOMATIC.M7
    rs6696489  6.592308e-21 2.743271e-76 2.002534e-167 3.562023e-275          0
    rs59788818 3.491208e-34 1.711087e-93 7.197813e-188 2.275374e-298          0
    rs2843493  1.037518e-11 5.167773e-67 3.034251e-158 3.725864e-266          0
                GERMLINE.M0  GERMLINE.M1  GERMLINE.M2  GERMLINE.M3   GERMLINE.M4
    rs6696489  1.723074e-04 9.989057e-01 6.216658e-06 3.302211e-27  1.067221e-82
    rs59788818 5.836546e-23 9.999757e-01 2.683928e-07 1.233256e-51 2.340805e-111
    rs2843493  8.802336e-49 1.023147e-07 9.997285e-01 2.798399e-18  3.718923e-74
                 GERMLINE.M5   GERMLINE.M6 GERMLINE.M7 GERMLINE.CONTHIGH
    rs6696489  6.678185e-174 1.070845e-281           0      3.060961e-62
    rs59788818 5.523414e-206 1.183050e-316           0     8.253575e-127
    rs2843493  9.760218e-166 6.969108e-274           0      5.841728e-59
               GERMLINE.CONTLOW ML.SOMATIC ML.M ML.C     ML.AR        AR
    rs6696489      2.105422e-22      FALSE    1    3 0.3571429 0.2803790
    rs59788818     3.106132e-86      FALSE    1    3 0.3571429 0.3957373
    rs2843493     4.644636e-138      FALSE    2    3 0.6428571 0.6259540
               CN.Subclonal Log.Ratio Prior.Somatic Prior.Contamination ML.LOH
    rs6696489         FALSE 1.5984124  0.0004950495                0.01  FALSE
    rs59788818        FALSE 0.2600908  0.0004950495                0.01  FALSE
    rs2843493         FALSE 0.7559721  0.0004950495                0.01  FALSE
               num.snps.loh.segment percentage.loh.in.loh.segment
    rs6696489                    49                             0
    rs59788818                   49                             0
    rs2843493                    49                             0

This lists all posterior probabilities for all possible SNV states.

