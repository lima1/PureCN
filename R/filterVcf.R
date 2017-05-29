#' Basic VCF filter function
#' 
#' Function to remove artifacts and low confidence/quality variant calls.
#' 
#' 
#' @param vcf \code{CollapsedVCF} object, read in with the \code{readVcf}
#' function from the VariantAnnotation package.
#' @param tumor.id.in.vcf The tumor id in the \code{CollapsedVCF} (optional).
#' @param use.somatic.status If somatic status and germline data is available,
#' then use this information to remove non-heterozygous germline SNPs or
#' germline SNPS with biased allelic fractions.
#' @param snp.blacklist A file with blacklisted genomic regions. Must
#' be parsable by \code{import} from \code{rtracklayer}, for a example a
#' BED file with file extension \sQuote{.bed}.
#' @param af.range Exclude SNPs with allelic fraction smaller or greater than
#' the two values, respectively. The higher value removes homozygous SNPs,
#' which potentially have allelic fractions smaller than 1 due to artifacts or
#' contamination. If a matched normal is available, this value is ignored,
#' because homozygosity can be confirmed in the normal.
#' @param contamination.range Count SNPs in dbSNP with allelic fraction
#' in the specified range. If the number of these putative contamination 
#' SNPs exceeds an expected value and if they are found on almost all 
#' chromosomes, the sample is flagged as potentially contaminated and extra
#' contamination estimation steps will be performed later on.
#' @param min.coverage Minimum coverage in tumor. Variants with lower coverage
#' are ignored.
#' @param min.base.quality Minimim base quality in tumor. Requires a \code{BQ}
#' genotype field in the VCF.
#' @param min.supporting.reads Minimum number of reads supporting the alt
#' allele.  If \code{NULL}, calculate based on coverage and assuming sequencing
#' error of 10^-3.
#' @param error Estimated sequencing error rate. Used to calculate minimum
#' number of supporting reads using \code{\link{calculatePowerDetectSomatic}}.
#' @param target.granges \code{GenomicRanges} object specifiying the target
#' postions. Used to remove off-target reads. If \code{NULL}, do not check
#' whether variants are on or off-target.
#' @param remove.off.target.snvs If set to a true value, will remove all SNVs
#' outside the covered regions.
#' @param model.homozygous If set to \code{TRUE}, does not remove homozygous
#' SNPs. Ignored in case a matched normal is provided in the VCF.
#' @param interval.padding Include variants in the interval flanking regions of
#' the specified size in bp. Requires \code{target.granges}.
#' @return A list with elements \item{vcf}{The filtered \code{CollapsedVCF}
#' object.} \item{flag}{A flag (\code{logical(1)}) if problems were
#' identified.} \item{flag_comment}{A comment describing the flagging.}
#' @author Markus Riester
#' @seealso \code{\link{calculatePowerDetectSomatic}}
#' @examples
#' 
#' # This function is typically only called by runAbsolute via 
#' # fun.filterVcf and args.filterVcf.
#' vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
#' vcf <- readVcf(vcf.file, "hg19")
#' vcf.filtered <- filterVcfBasic(vcf)        
#' 
#' @export filterVcfBasic
#' @importFrom GenomeInfoDb seqnames seqlevelsStyle seqlevelsStyle<-
#'             genomeStyles sortSeqlevels
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom stats pbeta
filterVcfBasic <- function(vcf, tumor.id.in.vcf = NULL, 
use.somatic.status = TRUE, snp.blacklist = NULL, af.range = c(0.03, 0.97),
contamination.range = c(0.01, 0.075), min.coverage = 15, min.base.quality = 25,
min.supporting.reads = NULL, error = 0.001, target.granges = NULL,
remove.off.target.snvs = TRUE, model.homozygous = FALSE, 
interval.padding = 50) {
    flag <- NA
    flag_comment <- NA

    if (is.null(tumor.id.in.vcf)) {
        tumor.id.in.vcf <- .getTumorIdInVcf(vcf)
    }
    if (use.somatic.status) {
        n <- nrow(vcf)
        vcf <- .testGermline(vcf, tumor.id.in.vcf)
        flog.info("Removing %i non heterozygous (in matched normal) germline SNPs.", 
            n-nrow(vcf))
    } else {
        info(vcf)$SOMATIC <- NULL
    }    

    # find supporting read cutoffs based on coverage and sequencing error
    #--------------------------------------------------------------------------
    if (is.null(min.supporting.reads)) {
        depths <- geno(vcf)$DP[,tumor.id.in.vcf]
        minDepth <- round(log2(mean(depths)))

        depths <- sort(unique(round(log2(depths))))
        depths <- depths[depths>=minDepth]
        depths <- 2^depths

        cutoffs <- sapply(depths, function(d) 
            calculatePowerDetectSomatic(d, ploidy=2, purity=1, error=error, 
                verbose=FALSE)$k)
        depths <- c(0,depths)
    } else { 
        depths <- 0
        cutoffs <- min.supporting.reads
    }
    n <- nrow(vcf)
    
    .sufficientReads <- function(vcf, ref, depths, cutoffs) {
        idx <- rep(TRUE, nrow(vcf))

        .filterVcfByAD <- function(vcf, min.supporting.reads, depth) {
            # remove variants with insufficient reads. This includes reference alleles
            # to remove homozygous germline variants.
            do.call(rbind, geno(vcf)$AD[,tumor.id.in.vcf])[,ifelse(ref,1,2)] >= 
                min.supporting.reads | 
                geno(vcf)$DP[,tumor.id.in.vcf] < depth
           }

        for (i in seq_along(cutoffs)) {
            idx <- idx & .filterVcfByAD(vcf, cutoffs[i], depths[i])
        }
        idx
    }    
    idxNotAlt <- .sufficientReads(vcf,ref=FALSE, depths, cutoffs)
    vcf <- vcf[idxNotAlt]
    idxNotHomozygous <- .sufficientReads(vcf,ref=TRUE, depths, cutoffs)

    #--------------------------------------------------------------------------

    # heurestic to find potential contamination
    #--------------------------------------------------------------------------
    idx <- info(vcf)$DB & idxNotHomozygous &
        ( unlist(geno(vcf)$FA[,tumor.id.in.vcf]) <  contamination.range[2] | 
        unlist(geno(vcf)$FA[,tumor.id.in.vcf]) > (1 - contamination.range[2]))

    fractionContaminated <- sum(idx)/sum(info(vcf)$DB & idxNotHomozygous)
    minFractionContaminated <- 0.02

    if (fractionContaminated > 0) {
        expectedAllelicFraction <- contamination.range[1] * 0.75
        powerDetectCont <- calculatePowerDetectSomatic(mean(geno(vcf)$DP[,tumor.id.in.vcf], 
            na.rm=TRUE), f=expectedAllelicFraction, verbose=FALSE)$power
        minFractionContaminated <- min(0.2, max(minFractionContaminated, powerDetectCont * 0.5))
        flog.info("Initial testing for significant sample cross-contamination: %s", 
            ifelse(fractionContaminated>minFractionContaminated, "maybe", "unlikely"))
    }

    # do we have many low allelic fraction calls that are in dbSNP on basically
    # all chromosomes? then we found some contamination
    if (sum(runLength(seqnames(rowRanges(vcf[idx])))>3) >= 20 &&
        fractionContaminated >= minFractionContaminated) {
        flag <- TRUE
        flag_comment <- "POTENTIAL SAMPLE CONTAMINATION"    
    }

    # If we have a matched normal, we can distinguish LOH from homozygous
    # in 100% pure samples. If not we need to see a sufficient number
    # of alt alleles to believe a heterozygous normal genotype.    
    if (use.somatic.status || model.homozygous) {
        af.range[2] <- 1
        if (model.homozygous) af.range[2] <- Inf
    } else {
        vcf <- vcf[idxNotHomozygous]
    }

    vcf <- vcf[unlist(geno(vcf)$DP[,tumor.id.in.vcf]) >= min.coverage]
    vcf <- vcf[unlist(geno(vcf)$FA[,tumor.id.in.vcf]) >= af.range[1]]
    # remove homozygous germline
    vcf <- vcf[!info(vcf)$DB | geno(vcf)$FA[,tumor.id.in.vcf] < af.range[2]]
    flog.info("Removing %i SNPs with AF < %.3f or AF >= %.3f or less than %i supporting reads or depth < %i.", 
        n-nrow(vcf), af.range[1], af.range[2], cutoffs[1], min.coverage)

    if (!is.null(snp.blacklist)) {
        for (i in seq_along(snp.blacklist)) {
            blackBed <- import(snp.blacklist[i])
            ov <- suppressWarnings(overlapsAny(vcf, blackBed))
            vcf <- vcf[!ov]
            flog.info("Removing %i blacklisted SNPs.", n-nrow(vcf))
        }    
    }

    if (!is.null(geno(vcf)$BQ)) {
        n.vcf.before.filter <- nrow(vcf)
        vcf <- vcf[which(as.numeric(geno(vcf)$BQ[,tumor.id.in.vcf])>=min.base.quality)]
        flog.info("Removing %i low quality variants with BQ < %i.", 
            n.vcf.before.filter - nrow(vcf), min.base.quality) 
    } else if (!is.null(rowRanges(vcf)$QUAL)) {
        n.vcf.before.filter <- nrow(vcf)
        idx <- which(as.numeric(rowRanges(vcf)$QUAL)>=min.base.quality)
        if (length(idx)) vcf <- vcf[idx]
        flog.info("Removing %i low quality variants with BQ < %i.", 
            n.vcf.before.filter - nrow(vcf), min.base.quality) 
    }     
    
    if (!is.null(target.granges)) {
        vcf <- .annotateVcfTarget(vcf, target.granges, interval.padding)
        if (remove.off.target.snvs) {
            n.vcf.before.filter <- nrow(vcf)
            # make sure all SNVs are in covered exons
            vcf <- vcf[info(vcf)$OnTarget>0]
            flog.info("Removing %i variants outside intervals.", 
                n.vcf.before.filter - nrow(vcf))
        }        
    }
    if (!is.null(info(vcf)$DB) && sum(info(vcf)$DB) < nrow(vcf)/2) {
        flog.warn("Less than half of variants in dbSNP. Make sure that VCF %s", 
            "contains both germline and somatic variants.")
    }
    list(
        vcf=vcf, 
        flag=flag, 
        flag_comment=flag_comment 
    )
}


#' Filter VCF MuTect
#' 
#' Function to remove artifacts and low confidence/quality calls from a MuTect
#' generated VCF file. Also applies filters defined in \code{filterVcfBasic}.
#' This function will only keep variants listed in the stats file and those not
#' matching the specified failure reasons.
#' 
#' 
#' @param vcf \code{CollapsedVCF} object, read in with the \code{readVcf}
#' function from the VariantAnnotation package.
#' @param tumor.id.in.vcf The tumor id in the VCF file, optional.
#' @param stats.file MuTect stats file.
#' @param ignore MuTect flags that mark variants for exclusion.
#' @param \dots Additional arguments passed to \code{\link{filterVcfBasic}}.
#' @return A list with elements \code{vcf}, \code{flag} and
#' \code{flag_comment}.  \code{vcf} contains the filtered \code{CollapsedVCF},
#' \code{flag} a \code{logical(1)} flag if problems were identified, further
#' described in \code{flag_comment}.
#' @author Markus Riester
#' @seealso \code{\link{filterVcfBasic}}
#' @examples
#' 
#' ### This function is typically only called by runAbsolute via the 
#' ### fun.filterVcf and args.filterVcf comments.
#' library(VariantAnnotation)    
#' vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
#' vcf <- readVcf(vcf.file, "hg19")
#' vcf.filtered <- filterVcfMuTect(vcf)        
#' 
#' @export filterVcfMuTect
filterVcfMuTect <- function(vcf, tumor.id.in.vcf = NULL, stats.file = NULL, 
ignore=c("clustered_read_position", "fstar_tumor_lod", "nearby_gap_events", 
"poor_mapping_region_alternate_allele_mapq", "poor_mapping_region_mapq0", 
"possible_contamination", "strand_artifact", "seen_in_panel_of_normals"),
... ){
    if (is.null(stats.file)) return(
        filterVcfBasic(vcf, tumor.id.in.vcf, ...))
    
    stats <- read.delim(stats.file, as.is=TRUE, skip=1)

    if (is.null(stats$contig) || is.null(stats$position)) {
        flog.warn("MuTect stats file lacks contig and position columns.")
        return(filterVcfBasic(vcf, tumor.id.in.vcf, ...))
    }    

    gr.stats <- GRanges(seqnames=stats$contig, 
        IRanges(start=stats$position, end=stats$position))
    
    ov <- findOverlaps(vcf, gr.stats)
     
    if (!identical(queryHits(ov),subjectHits(ov)) || 
            nrow(vcf) != nrow(stats)) {
        n <- nrow(vcf)
        stats <- stats[subjectHits(ov),]
        vcf <- vcf[queryHits(ov)]
        flog.warn("MuTect stats file and VCF file do not align perfectly. Will remove %i unmatched variants.",
            n-nrow(vcf))
    }    
    if (is.null(stats$failure_reasons)) {
        flog.warn("MuTect stats file lacks failure_reasons column.%s",
            " Keeping all variants listed in stats file.")
        return(filterVcfBasic(vcf, tumor.id.in.vcf, ...))
    }

    n <- nrow(vcf)

    ids <- sort(unique(unlist(sapply(ignore, grep, stats$failure_reasons))))
    vcf <- vcf[-ids]

    flog.info("Removing %i MuTect calls due to blacklisted failure reasons.", 
        n-nrow(vcf))
    filterVcfBasic(vcf, tumor.id.in.vcf, ...)
}


#' Set Somatic Prior VCF
#' 
#' Function to set prior for somatic mutation status for each variant in the
#' provided \code{CollapsedVCF} object.
#' 
#' 
#' @param vcf \code{CollapsedVCF} object, read in with the \code{readVcf}
#' function from the VariantAnnotation package.
#' @param prior.somatic Prior probabilities for somatic mutations. First value
#' is for the case when no matched normals are available and the variant is not
#' in dbSNP (second value). Third value is for variants with MuTect somatic
#' call. Different from 1, because somatic mutations in segments of copy number
#' 0 have 0 probability and artifacts can thus have dramatic influence on
#' likelihood score. Forth value is for variants not labeled as somatic by
#' MuTect. Last two values are optional, if vcf contains a flag Cosmic.CNT, it
#' will set the prior probability for variants with CNT > 2 to the first of
#' those values in case of no matched normal available (0.995 default).  Final
#' value is for the case that variant is in both dbSNP and COSMIC > 2.
#' @param tumor.id.in.vcf Id of tumor in case multiple samples are stored in
#' VCF.
#' @param min.cosmic.cnt Minimum number of hits in the COSMIC database to 
#' call variant as likely somatic.
#' @return A \code{numeric(nrow(vcf))} vector with the prior probability of
#' somatic status for each variant in the \code{CollapsedVCF}.
#' @author Markus Riester
#' @examples
#' 
#' # This function is typically only called by runAbsoluteCN via the 
#' # fun.setPriorVcf and args.setPriorVcf comments.
#' vcf.file <- system.file("extdata", "example_vcf.vcf", package="PureCN")
#' vcf <- readVcf(vcf.file, "hg19")
#' vcf.priorsomatic <- setPriorVcf(vcf)        
#' 
#' @export setPriorVcf
setPriorVcf <- function(vcf,
prior.somatic=c(0.5, 0.0005, 0.999, 0.0001, 0.995, 0.5), 
tumor.id.in.vcf=NULL,
min.cosmic.cnt=4) {
    if (is.null(tumor.id.in.vcf)) {
        tumor.id.in.vcf <- .getTumorIdInVcf(vcf)
    }
    if (!is.null(info(vcf)$SOMATIC)) {
         tmp <- prior.somatic
         prior.somatic <- ifelse(info(vcf)$SOMATIC,
            prior.somatic[3],prior.somatic[4])

         flog.info("Found SOMATIC annotation in VCF.")
         flog.info("Setting somatic prior probabilities for somatic variants to %f or to %f otherwise.", 
            tmp[3], tmp[4])
    } else {
         tmp <- prior.somatic
         prior.somatic <- ifelse(info(vcf)$DB,
            prior.somatic[2], prior.somatic[1])
         if (!is.null(info(vcf)$Cosmic.CNT)) {
             flog.info("Found COSMIC annotation in VCF.")
             flog.info("Setting somatic prior probabilities for hits to %f or to %f if in both COSMIC and dbSNP.", 
                tmp[5], tmp[6])

             prior.somatic[which(info(vcf)$Cosmic.CNT >= min.cosmic.cnt)] <- tmp[5]
             prior.somatic[which(info(vcf)$Cosmic.CNT >= min.cosmic.cnt & 
                info(vcf)$DB)] <- tmp[6]
         } else {
             flog.info("Setting somatic prior probabilities for dbSNP hits to %f or to %f otherwise.", 
                tmp[2], tmp[1])
         }      
    }     
    prior.somatic
}

# If a matched normal is available, this will remove all
# heterozygous germline SNPs with biased allelic fractions   
.testGermline <-
function(vcf, tumor.id.in.vcf, allowed=0.05) {
    # extract normal allelic fractions and total coverage from the VCF
    arAll <- as.numeric(geno(vcf)$FA[,-match(tumor.id.in.vcf,
        colnames(geno(vcf)$FA))])
    dpAll <- as.numeric(geno(vcf)$DP[,-match(tumor.id.in.vcf, 
        colnames(geno(vcf)$DP))])
    
    # calculate probability that true allelic fraction in normal is 0.5
    pBeta <-pbeta(1/2, 
            shape1=arAll*dpAll+1, 
            shape2=(1-arAll)*dpAll+1,log.p=FALSE)
    
    # keep only somatic, non-biased and if allelic ratio is close
    # enough. The latter is useful for ultra-deep sequencing, when
    # non-reference bias can lead to small p-values.
    idx <- info(vcf)$SOMATIC | (pBeta > 0.025 & pBeta < 1-0.025) | 
           (arAll > 0.5-allowed & arAll < 0.5+allowed) 
    idx <- idx & dpAll > 5       
    vcf[idx]
}
# calculate allelic fraction from read depths
.addFaField <- function(vcf, field="FA") {
    if (is.null(geno(vcf)$AD)) {
        .stopRuntimeError("No AD geno field in VCF.")
    }
    matrixFA <- do.call(rbind, apply(geno(vcf)$AD,1, function(x) lapply(x,
        function(y) y[2]/sum(y))))
    newGeno <- DataFrame(
        Number="A", Type="Float",
        Description="Allele fraction of the alternate allele",
        row.names=field)
    geno(header(vcf)) <- rbind(geno(header(vcf)), newGeno)
    geno(vcf)[[field]] <- matrixFA
    vcf
}
.addDpField <- function(vcf, field="DP") {
    if (is.null(geno(vcf)$AD)) {
        .stopRuntimeError("No AD geno field in VCF.")
    }
    matrixDP <- apply(geno(vcf)$AD,2,sapply, sum)
    newGeno <- DataFrame(
        Number=1, Type="Integer",
        Description="Approximate read depth",
        row.names=field)
    geno(header(vcf)) <- rbind(geno(header(vcf)), newGeno)
    geno(vcf)[[field]] <- matrixDP
    vcf
}
.getNormalIdInVcf <- function(vcf, tumor.id.in.vcf) {
    samples(header(vcf))[-match(tumor.id.in.vcf, samples(header(vcf)))] 
}    
.getTumorIdInVcf <- function(vcf, sampleid=NULL) {
    .getTumorId <- function(vcf, sampleid) {
        if (ncol(vcf) == 1) return(1)
        if (ncol(vcf) != 2) {
            .stopUserError("VCF contains ", ncol(vcf), " samples. Should be ",
                "either one or two.")
        }        
        if (!is.null(geno(vcf)$GT)) {
            cs <- colSums(geno(vcf)$GT=="0")
            if (max(cs) > 0) return(which.min(cs))
            cs <- colSums(geno(vcf)$GT=="0/0")
            if (max(cs) > 0) return(which.min(cs))
        }
        if (!is.null(geno(vcf)$FA)) {
            cs <- apply(geno(vcf)$FA,2,function(x) sum(as.numeric(x)==0))
            if (max(cs) > 0) return(which.min(cs))
        }    
        if (sampleid %in% samples(header(vcf))) {
            return(match(sampleid, samples(header(vcf))))
        }    
        flog.warn("Cannot determine tumor vs. normal in VCF. %s",
            "Assuming it is first sample. Specify tumor via sampleid argument.")
        return(1)
    }
    tumor.id.in.vcf <- .getTumorId(vcf, sampleid=sampleid)
    if (!is.numeric(tumor.id.in.vcf) || tumor.id.in.vcf < 1 ||
         tumor.id.in.vcf > 2) {
        .stopRuntimeError("Tumor id not in expected range.")
    }    
    tumor.id.in.vcf <- samples(header(vcf))[tumor.id.in.vcf]    
    # check that sampleid matches tumor and not normal.
    if (!is.null(sampleid)) {
        if (sampleid %in% samples(header(vcf)) && 
            sampleid != tumor.id.in.vcf) {
            flog.warn("Sampleid looks like a normal in VCF, not like a tumor.")
        }    
    }    
    tumor.id.in.vcf
}
.checkVcfFieldAvailable <- function(vcf, field) {
    if (is.null(geno(vcf)[[field]]) ||  
        sum(is.finite(as.numeric(geno(vcf)[[field]][,1])))<1) {
        return(FALSE)
    }    
    return(TRUE)
}    
.readAndCheckVcf <- function(vcf.file, genome) {
    if (class(vcf.file) == "character") {    
        vcf <- readVcf(vcf.file, genome)
    } else if (class(vcf.file) != "CollapsedVCF") {
        .stopUserError("vcf.file neither a filename nor a CollapsedVCF ", 
            "object.") 
    } else {
        vcf <- vcf.file
    } 
    triAllelic <- elementNROWS(alt(vcf))>1
    if (sum(triAllelic)) {
        n <- nrow(vcf)
        vcf <- vcf[which(!triAllelic)]
        flog.info("Removing %i triallelic sites.",n-nrow(vcf))
    }    
    if (is.null(info(vcf)$DB)) {
        # try to add an DB field based on rownames
        vcf <- .addDbField(vcf)
    }
    if (is.null(geno(vcf)$AD)) {
        vcf <- .addADField(vcf)
    }
    if (!.checkVcfFieldAvailable(vcf, "FA")) {
        # try to add an FA geno field if missing
        vcf <- .addFaField(vcf)
    }
    if (!.checkVcfFieldAvailable(vcf, "DP")) {
        # try to add an DP geno field if missing
        vcf <- .addDpField(vcf)
    }
    vcf     
}    
.addDbField <- function(vcf) {
     db <- grepl("^rs",rownames(vcf))
     if (!sum(db)) {
        .stopUserError("vcf.file has no DB info field for dbSNP membership.")
     } else  { 
        flog.warn("vcf.file has no DB info field for dbSNP membership.%s",
            " Guessing it based on ID.")
     }   
    newInfo <- DataFrame(
        Number=0, Type="Flag",
        Description="dbSNP Membership",
        row.names="DB")
    info(header(vcf)) <- rbind(info(header(vcf)), newInfo)
    info(vcf)$DB <- db
    vcf
}

.addADField <- function(vcf, field="AD") {
    # FreeBayes
    if (!length(setdiff(c("DP","AO", "RO"),names(geno(vcf))))) {
        matrixAD <- do.call(cbind, lapply(samples(header(vcf)), function(j) { 
            ao <- unlist(geno(vcf)$AO[,j])
            ro <- unlist(geno(vcf)$RO[,j])
            AD <- lapply(seq_along(ao), function(i) as.integer(c(ro[i], ao[i])))
            names(AD) <- names(ao)
            AD
        })) 
        colnames(matrixAD) <- samples(header(vcf))   
        newGeno <- DataFrame(
            Number=".", Type="Integer",
            Description="Allelic depths for the ref and alt alleles in the order listed",
            row.names=field)
        geno(header(vcf)) <- rbind(geno(header(vcf)), newGeno)
        geno(vcf)[[field]] <- matrixAD
    } else {   
        .stopUserError("vcf.file has no AD geno field containing read depths ",
            "of ref and alt.")
    }
    vcf
}
    
.addCosmicCNT <- function(vcf, cosmic.vcf.file) {
    if (!is.null(info(vcf)$Cosmic.CNT)) {
        flog.info("VCF already COSMIC annotated. Skipping.")
        return(vcf)        
    }
    cosmicSeqStyle <- seqlevelsStyle(headerTabix(
        TabixFile(cosmic.vcf.file))$seqnames)

    vcfRenamedSL <- vcf
    if (!length(intersect(seqlevelsStyle(vcf), cosmicSeqStyle))) {
        seqlevelsStyle(vcfRenamedSL) <- cosmicSeqStyle[1]
    }
    flog.info("Reading COSMIC VCF...")
    # look-up the variants in COSMIC
    cosmic.vcf <- readVcf(cosmic.vcf.file, genome=genome(vcf)[1],  
        ScanVcfParam(which = rowRanges(vcfRenamedSL),
            info="CNT",
            fixed="ALT",
            geno=NA
        )
    )
    ov <- findOverlaps(vcfRenamedSL, cosmic.vcf, type="equal")
    # make sure that alt alleles match
    idx <- as.logical(alt(vcf[queryHits(ov)]) ==
        alt(cosmic.vcf[subjectHits(ov)]))
    ov <- ov[idx]
    if (!length(ov)) return(vcf)

    if (is.null(info(cosmic.vcf)$CNT)) {
        flog.warn("Cosmic VCF has no CNT info field. %s",
            "Giving up COSMIC annotation.")
        return(vcf)
    }
    
    newInfo <- DataFrame(
        Number=1, Type="Integer",
        Description="How many samples in Cosmic have this mutation",
        row.names="Cosmic.CNT")
    info(header(vcf)) <- rbind(info(header(vcf)), newInfo)
    info(vcf)$Cosmic.CNT <- NA
    info(vcf)$Cosmic.CNT[queryHits(ov)] <- info(cosmic.vcf[subjectHits(ov)])$CNT
    vcf
}

.calcTargetedGenome <- function(granges) {
    tmp <- reduce(granges)
    sum(width(tmp))/(1000^2)
}    

.annotateVcfTarget <- function(vcf, target.granges, interval.padding) {
    target.granges.padding <- .padGranges(target.granges, interval.padding)

    flog.info("Total size of targeted genomic region: %.2fMb (%.2fMb with %ibp padding).", 
        .calcTargetedGenome(target.granges), 
        .calcTargetedGenome(target.granges.padding),
        interval.padding)

    idxTarget <- overlapsAny(vcf, target.granges)
    idxPadding <- overlapsAny(vcf, target.granges.padding)
   
    newInfo <- DataFrame(
        Number=1, Type="Integer",
        Description="1: On-target; 2: Flanking region; 0: Off-target.",
        row.names="OnTarget")
    
    info(header(vcf)) <- rbind(info(header(vcf)), newInfo)
    info(vcf)$OnTarget <- 0
    info(vcf)$OnTarget[idxPadding] <- 2
    info(vcf)$OnTarget[idxTarget] <- 1
    
    # report stats in log file
    targetsWithSNVs <- overlapsAny(target.granges.padding, vcf)
    percentTargetsWithSNVs <- sum(targetsWithSNVs,na.rm=TRUE)/
        length(targetsWithSNVs)*100
    tmp <- ""
    if (percentTargetsWithSNVs > 20) { 
        tmp <- " segmentationPSCBS might produce better results."
    }
    flog.info("%.1f%% of targets contain SNVs. %s", 
        percentTargetsWithSNVs, tmp)

    vcf
}
