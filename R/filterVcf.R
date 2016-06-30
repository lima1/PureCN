filterVcfBasic <-
structure(function(#Basic VCF filter function
### Function to remove artifacts and low confidence/quality 
### variant calls.
vcf, 
### CollapsedVCF object, read in with the readVcf function 
### from the VariantAnnotation package.
tumor.id.in.vcf=NULL, 
### The tumor id in the CollapsedVCF (optional).
use.somatic.status=TRUE, 
### If somatic status and germline data is available, then use this 
### information to remove non-heterozygous germline SNPs or germline SNPS 
### with biased allelic fractions.
snp.blacklist=NULL,
### CSV file with SNP ids with expected allelic fraction 
### significantly different from 0.5 in diploid genomes. Can be an array of 
### lists. The function createSNPBlacklist can provide appropriate black 
### lists. Can also be a BED file (either tab or comma separated) of 
### blacklisted genomic regions (columns 1-3: chromosome, start, end).
af.range=c(0.03, 0.97),
### Exclude SNPs with allelic fraction smaller or greater than the 
### two values, respectively.
contamination.cutoff=c(0.05,0.075),
### Count SNPs in dbSNP with allelic fraction smaller than the 
### first value, if found on most chromosomes, remove all with AF smaller than
### the second value.
coverage.cutoff=20,
### Minimum coverage in tumor. Variants with lower coverage are ignored.
min.supporting.reads=NULL,
### Minimum number of reads supporting the alt allele. 
### If NULL, calculate based on coverage and assuming sequencing error 
### of 10^-3.
error=0.001,
### Estimated sequencing error rate. Used to calculate minimum
### number of supporting reads using calculatePowerDetectSomatic.
verbose=TRUE
) {
    flag <- NA
    flag_comment <- NA

    if (is.null(tumor.id.in.vcf)) {
        tumor.id.in.vcf <- names( which.min(colSums(geno(vcf)$GT=="0")) )
    }
    if (use.somatic.status) {
        n <- nrow(vcf)
        vcf <- .testGermline(vcf, tumor.id.in.vcf)
        if (verbose) message(paste("Removing", n-nrow(vcf), 
            "non heterozygous (in matched normal) germline SNPs."))
    } else {
        info(vcf)$SOMATIC <- NULL
    }    

    if (is.null(min.supporting.reads)) {
        min.supporting.reads <- calculatePowerDetectSomatic(
            mean(geno(vcf)$DP[,tumor.id.in.vcf], na.rm=TRUE),
            purity=1,ploidy=2, error=error, verbose=FALSE)$k
    }    

    n <- nrow(vcf)
    
    # remove variants with insufficient reads. This includes reference alleles
    # to remove homozygous germline variants.
    vcf <- vcf[do.call(rbind, 
        geno(vcf)$AD[,tumor.id.in.vcf])[,2] >= min.supporting.reads]
    vcf <- vcf[do.call(rbind, 
        geno(vcf)$AD[,tumor.id.in.vcf])[,1] >= min.supporting.reads]

    vcf <- vcf[unlist(geno(vcf)$DP[,tumor.id.in.vcf]) >= coverage.cutoff]
    vcf <- vcf[unlist(geno(vcf)$FA[,tumor.id.in.vcf]) >= af.range[1]]
    # remove homozygous germline
    vcf <- vcf[!info(vcf)$DB | geno(vcf)$FA[,tumor.id.in.vcf] < af.range[2]]
    if (verbose) message("Removing ", n-nrow(vcf), 
        " SNPs with AF < ", af.range[1],
        " or AF >= ", af.range[2], 
        " or less than ", min.supporting.reads, 
        " supporting reads or depth < ",
        coverage.cutoff, ".")

    if (!is.null(snp.blacklist)) {
        for (i in 1:length(snp.blacklist)) {
            snp.blacklist.data <- read.csv(snp.blacklist[i], as.is=TRUE)
            snp.blacklist.data2 <- read.delim(snp.blacklist[i], as.is=TRUE)
            if (ncol(snp.blacklist.data2) > ncol(snp.blacklist.data)) {
                snp.blacklist.data <- snp.blacklist.data2
            }
            n <- nrow(vcf)
            if (sum( rownames(vcf) %in% snp.blacklist.data[,1]) > 1  ) {
                vcf <- vcf[!rownames(vcf) %in% snp.blacklist.data[,1],]
            } else {
                if (is.null(snp.blacklist.data$seg.mean)) {
                    snp.blacklist.data$seg.mean <- 1
                }    
                snp.blacklist.data <- 
                    snp.blacklist.data[snp.blacklist.data$seg.mean > 0.2,]
                ov <- suppressWarnings(findOverlaps(vcf, 
                    GRanges(seqnames=snp.blacklist.data[,1], 
                    IRanges(start=snp.blacklist.data[,2], 
                            end=snp.blacklist.data[,3]))))

                idx <- !(1:nrow(vcf) %in% queryHits(ov) & info(vcf)$DB)
                vcf <- vcf[idx]
            }    
            if (verbose) message("Removing ", n-nrow(vcf), 
                " blacklisted SNPs.")
        }    
    }    

    idx <- info(vcf)$DB & unlist(geno(vcf)$FA[,tumor.id.in.vcf]) < 
        contamination.cutoff[1]

    # do we have many low allelic fraction calls that are in dbSNP on basically
    # all chromosomes? then we found some contamination
    if (sum(runLength(seqnames(rowRanges(vcf[idx])))>1) >= 20) {
        idx <- info(vcf)$DB & unlist(geno(vcf)$FA[,tumor.id.in.vcf]) < 
            contamination.cutoff[2]
        vcf <- vcf[which(!idx)]
        if (verbose) message("Removing ", sum(idx, na.rm=TRUE), 
            " contamination SNPs.")
        flag <- TRUE
        flag_comment <- "POTENTIAL SAMPLE CONTAMINATION"    
    }
    
    # if we have a matched normal, we can filter homozygous germline SNPs
    # easily, if not, we still should remove most of them by removing everything
    # with AF >= 0.97. But sometimes samples are very noisy, for example due to
    # contamination, and homozygous germline SNPs are between 0.9 and 1.  This
    # code finds samples which have a lot of dbSNPs between 0.95 and 0.97 on
    # almost all chromosomes.  If found, then remove everything that's in dbSNP
    # above 0.9.
    if (!use.somatic.status) {    
        idx <- info(vcf)$DB & unlist(geno(vcf)$FA[,tumor.id.in.vcf]) > 
            (1 - contamination.cutoff[1])

        # do we have many high allelic fraction calls that are in dbSNP on
        # basically all chromosomes?  then we found a very noisy sample
        if (sum(runLength(seqnames(rowRanges(vcf[idx])))>2) >= 20) {
            idx <- info(vcf)$DB & unlist(geno(vcf)$FA[,tumor.id.in.vcf]) > 
                (1 - contamination.cutoff[2])
            vcf <- vcf[which(!idx)]
            if (verbose) message("Removing ", sum(idx, na.rm=TRUE), 
                " noisy homozygous germline SNPs.")
            flag <- TRUE
            flag_comment <- .appendComment(flag_comment, 
                "NOISY HOMOZYGOUS GERMLINE CALLS")    
        }    
    }

    ##value<< A list with elements
    list(
        vcf=vcf, ##<< The filtered CollapsedVCF object.
        flag=flag, ##<< A flag (TRUE/FALSE) if problems were identified.
        flag_comment=flag_comment ##<< A comment describing the flagging. 
    )
##end<<
}, ex=function() {
# This function is typically only called by runAbsolute via the 
# fun.filterVcf and args.filterVcf comments.
library(VariantAnnotation)    
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
vcf <- readVcf(vcf.file, "hg19")
vcf.filtered <- filterVcfBasic(vcf)        
})    

filterVcfMuTect <- structure(function(#Filter VCF MuTect
### Function to remove artifacts and low confidence/quality calls from
### a MuTect generated VCF file.
vcf, 
### VCF object, read in with the readVcf function from the 
### VariantAnnotation package.
tumor.id.in.vcf=NULL, 
### The tumor id in the VCF file, optional.
stats.file=NULL, 
### MuTect stats file
ignore=c("clustered_read_position", "fstar_tumor_lod", "nearby_gap_events", 
"poor_mapping_region_alternate_allele_mapq", "poor_mapping_region_mapq0", 
"possible_contamination", "strand_artifact"),
### Failure flags that lead to exclusion of variant.
verbose=TRUE,
### Verbose output.
...
### Additional arguments passed to filterVcfBasic
){
    if (is.null(stats.file)) return(
        filterVcfBasic(vcf, tumor.id.in.vcf, verbose=verbose, ...))
    
    stats <- read.delim(stats.file, as.is=TRUE, skip=1)
    gr.stats <- GRanges(seqnames=stats$contig, 
        IRanges(start=stats$position, end=stats$position))
    
    ov <- findOverlaps(vcf, gr.stats)
     
    if (!identical(queryHits(ov),subjectHits(ov))) {
        warning("MuTect stats file and VCF file do not align perfectly. ",
         "Will remove unmatched variants.")
        stats <- stats[subjectHits(ov),]
        vcf <- vcf[queryHits(ov)]
    }    

    n <- nrow(vcf)

    ids <- sort(unique(unlist(sapply(ignore, grep, stats$failure_reasons))))
    vcf <- vcf[-ids]

    if (verbose) message("Removing ", n-nrow(vcf), 
        " MuTect calls due to blacklisted failure reasons.")
    filterVcfBasic(vcf, tumor.id.in.vcf, verbose=verbose, ...)
### A list with elements vcf, flag and flag_comment. "vcf" contains the 
### filtered CollapsedVCF, "flag" a flag if problems were identified, further 
### described in "flag_comment".    
},ex=function() {
### This function is typically only called by runAbsolute via the 
### fun.filterVcf and args.filterVcf comments.
library(VariantAnnotation)    
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
vcf <- readVcf(vcf.file, "hg19")
vcf.filtered <- filterVcfMuTect(vcf)        
})

setPriorVcf <- structure(function(# Set Somatic Prior VCF
### Function to set prior for somatic mutation status for each
### variant in the provided CollapsedVCF object.
vcf,
### CollapsedVCF object, read in with the readVcf function 
### from the VariantAnnotation package.
prior.somatic=c(0.5, 0.0005, 0.999, 0.0001, 0.995, 0.01), 
### Prior probabilities for somatic mutations. First value is for 
### the case when no matched normals are available and the variant is not in 
### dbSNP (second value). Third value is for variants with MuTect somatic call.
### Different from 1, because somatic mutations in segments of copy number 0 
### have 0 probability and artifacts can thus have dramatic influence on 
### likelihood score. Forth value is for variants not labeled as somatic by 
### MuTect. Last two values are optional, if vcf contains a flag Cosmic.CNT, 
### it will set the prior probability for variants with CNT > 2 to the first 
### of those values in case of no matched normal available (0.995 default). 
### Final value is for the case that variant is in both dbSNP and COSMIC > 2. 
verbose=TRUE
### Verbose output.
) {
    if (!is.null(info(vcf)$SOMATIC)) {
         tmp <- prior.somatic
         prior.somatic <- ifelse(info(vcf)$SOMATIC,
            prior.somatic[3],prior.somatic[4])
         if (verbose) message("Found SOMATIC annotation in VCF. ",
            "Setting somatic prior probabilities for somatic variants to ", 
            tmp[3]," or to ", tmp[4], " otherwise.")
    } else {
         tmp <- prior.somatic
         prior.somatic <- ifelse(info(vcf)$DB,
            prior.somatic[2], prior.somatic[1])
         if (!is.null(info(vcf)$Cosmic.CNT)) {
             if (verbose) message("Found COSMIC annotation in VCF. ",
                "Setting somatic prior probabilities for hits to ", tmp[5],
                " or to ", tmp[6], " if in both COSMIC and dbSNP.")

             prior.somatic[which(info(vcf)$Cosmic.CNT>2)] <- tmp[5]
             prior.somatic[which(info(vcf)$Cosmic.CNT>2 & 
                info(vcf)$DB)] <- tmp[6]
         } else {
             if (verbose) message("Setting somatic prior probabilities ",
                "for dbSNP hits to ", tmp[2]," or to ", tmp[1], " otherwise.")
         }      
    }     
    prior.somatic
### A vector with the prior probability of somatic status for each 
### variant in the CollapsedVCF.
},ex=function() {
# This function is typically only called by runAbsoluteCN via the 
# fun.setPriorVcf and args.setPriorVcf comments.
vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
vcf <- readVcf(vcf.file, "hg19")
vcf.priorsomatic <- setPriorVcf(vcf)        
})     
   
.testGermline <-
function(vcf, tumor.id.in.vcf) {
    ar_all <- do.call(rbind, geno(vcf)$FA[,-match(tumor.id.in.vcf, 
        colnames(geno(vcf)$FA)), drop=FALSE])
    dp_all <- geno(vcf)$DP[,-match(tumor.id.in.vcf, 
        colnames(geno(vcf)$DP)), drop=FALSE]
    xx <-sapply(1:length(ar_all), function(j) 
        pbeta(1/2, 
            shape1=ar_all[j]*dp_all[j]+1, 
            shape2=(1-ar_all[j])*dp_all[j]+1,log.p=FALSE))
    idx <- info(vcf)$SOMATIC | (xx > 0.025 & xx < 1-0.025)
    vcf[idx]
}
