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
#' @param af.range Exclude variants with allelic fraction smaller or greater than
#' the two values, respectively. The higher value removes homozygous SNPs,
#' which potentially have allelic fractions smaller than 1 due to artifacts or
#' contamination. If a matched normal is available, this value is ignored,
#' because homozygosity can be confirmed in the normal.
#' @param contamination.range Count variants in dbSNP with allelic fraction
#' in the specified range. If the number of these putative contamination 
#' variants exceeds an expected value and if they are found on almost all 
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
#' variants. Ignored in case a matched normal is provided in the VCF.
#' @param interval.padding Include variants in the interval flanking regions of
#' the specified size in bp. Requires \code{target.granges}.
#' @param DB.info.flag Flag in INFO of VCF that marks presence in common
#' germline databases. Defaults to \code{DB} that may contain somatic variants
#' if it is from an unfiltered dbSNP VCF.
#' @return A list with elements \item{vcf}{The filtered \code{CollapsedVCF}
#' object.} \item{flag}{A flag (\code{logical(1)}) if problems were
#' identified.} \item{flag_comment}{A comment describing the flagging.}
#' @author Markus Riester
#' @seealso \code{\link{calculatePowerDetectSomatic}}
#' @examples
#' 
#' # This function is typically only called by runAbsolute via 
#' # fun.filterVcf and args.filterVcf.
#' vcf.file <- system.file("extdata", "example.vcf.gz", package="PureCN")
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
interval.padding = 50, DB.info.flag = "DB") {
    flag <- NA
    flag_comment <- NA

    if (is.null(tumor.id.in.vcf)) {
        tumor.id.in.vcf <- .getTumorIdInVcf(vcf)
    }
    if (use.somatic.status) {
        n <- .countVariants(vcf)
        vcf <- .testGermline(vcf, tumor.id.in.vcf)
        flog.info("Removing %i non heterozygous (in matched normal) germline SNPs.", 
            n - .countVariants(vcf))
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
    n <- .countVariants(vcf)
        
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
    vcf <- .removeVariants(vcf, !idxNotAlt, "sufficient reads")
    idxNotHomozygous <- .sufficientReads(vcf,ref=TRUE, depths, cutoffs)
    if (!sum(idxNotHomozygous)) .stopUserError("None of the heterozygous variants in provided VCF passed filtering.")
    #--------------------------------------------------------------------------

    # heurestic to find potential contamination
    #--------------------------------------------------------------------------
    idx <- info(vcf)[[DB.info.flag]] & idxNotHomozygous &
        (unlist(geno(vcf)$FA[,tumor.id.in.vcf]) < contamination.range[2] |
        unlist(geno(vcf)$FA[,tumor.id.in.vcf]) > (1 - contamination.range[2]))
    
    if (!sum(idx)) {
        # this usually only happens with wrong input where all germline SNPs are removed
        fractionContaminated <- 0
    } else {
        fractionContaminated <- sum(idx)/sum(info(vcf)[[DB.info.flag]] & idxNotHomozygous)
    }
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
        vcf <- .removeVariants(vcf, !idxNotHomozygous, "homozygous")
    }

    vcf <- .removeVariants(vcf, 
        unlist(geno(vcf)$DP[,tumor.id.in.vcf]) < min.coverage,
        "min.coverage")

    vcf <- .removeVariants(vcf, 
        unlist(geno(vcf)$FA[,tumor.id.in.vcf]) < af.range[1],
        "af.range")

    # remove homozygous germline
    vcf <- .removeVariants(vcf, 
        info(vcf)[[DB.info.flag]] &
        geno(vcf)$FA[,tumor.id.in.vcf] >= af.range[2],
        "homozygous af.range")

    flog.info("Removing %i variants with AF < %.3f or AF >= %.3f or less than %i supporting reads or depth < %i.", 
        n-.countVariants(vcf), af.range[1], af.range[2], cutoffs[1], min.coverage)
    n <- .countVariants(vcf)

    if (!is.null(snp.blacklist)) {
        for (i in seq_along(snp.blacklist)) {
            blackBed <- try(import(snp.blacklist[i]))
            if (class(blackBed) == "try-error") {
                .stopUserError("Could not import snp.blacklist ", snp.blacklist[[i]],
                    ":", blackBed)
            }    
            ov <- suppressWarnings(overlapsAny(vcf, blackBed))
            vcf <- .removeVariants(vcf, ov, "blacklist")
            flog.info("Removing %i blacklisted variants.", 
                      n-.countVariants(vcf))
        }    
    }
    
    if (!is.null(min.base.quality) && min.base.quality > 0) {
        vcf <- .filterVcfByBQ(vcf, tumor.id.in.vcf, min.base.quality)
    }

    if (!is.null(target.granges)) {
        vcf <- .annotateVcfTarget(vcf, target.granges, interval.padding)
        if (remove.off.target.snvs) {
            n.vcf.before.filter <- .countVariants(vcf)
            # make sure all SNVs are in covered exons
            key <- paste0(.getPureCNPrefixVcf(vcf), "OnTarget")
            vcf <- .removeVariants(vcf, info(vcf)[[key]] <= 0, "intervals") 
            flog.info("Removing %i variants outside intervals.", 
                n.vcf.before.filter - .countVariants(vcf))
        }
    }
    if (!is.null(info(vcf)[[DB.info.flag]]) && sum(info(vcf)[[DB.info.flag]]) < nrow(vcf)/2) {
        flog.warn("Less than half of variants in dbSNP. Make sure that VCF %s", 
            "contains both germline and somatic variants.")
    }
    list(
        vcf=vcf, 
        flag=flag, 
        flag_comment=flag_comment 
    )
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
            shape1 = arAll*dpAll+1,
            shape2 = (1-arAll)*dpAll+1,log.p=FALSE)
    
    # keep only somatic, non-biased and if allelic ratio is close
    # enough. The latter is useful for ultra-deep sequencing, when
    # non-reference bias can lead to small p-values.
    idx <- info(vcf)$SOMATIC | (pBeta > 0.025 & pBeta < 1-0.025) | 
           (arAll > 0.5-allowed & arAll < 0.5+allowed) 
    idx <- idx & dpAll > 5       
    .removeVariants(vcf, !idx, "matched germline")
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
.addFilterFlag <- function(vcf) {
    newFilter <- DataFrame(
        Description = "Ignored by PureCN",
        row.names = "purecn_ignore")
    fixed(header(vcf))$FILTER <- rbind(fixed(header(vcf))$FILTER, newFilter)
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
            cs <- apply(geno(vcf)$FA,2,function(x) length(which(as.numeric(x)==0)))
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
            tumor.id.in.vcf <- sampleid
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
.readAndCheckVcf <- function(vcf.file, genome, DB.info.flag = "DB", 
                             POPAF.info.field = "POP_AF", 
                             min.pop.af = 0.001, check.DB = TRUE,
                             vcf.field.prefix = NULL) {
    if (is(vcf.file, "character")) {
        vcf <- readVcf(vcf.file, genome)
    } else if (!is(vcf.file, "CollapsedVCF")) {
        .stopUserError("vcf.file neither a filename nor a CollapsedVCF ", 
            "object.") 
    } else {
        vcf <- vcf.file
    }
    # add PureCN ignore flag
    vcf <- .addFilterFlag(vcf) 
    flog.info("Found %i variants in VCF file.", length(vcf))
    triAllelic <- elementNROWS(alt(vcf))>1
    if (sum(triAllelic)) {
        n <- .countVariants(vcf)
        vcf <- .removeVariants(vcf, triAllelic, "triallelic")
        flog.info("Removing %i triallelic sites.",n-length(vcf))
    } 
    if (is.null(info(vcf)$SOMATIC)) {
        # try to add a SOMATIC flag for callers that do not
        # provide one (matched tumor/normal cases only)
        vcf <- .addSomaticField(vcf)
    }
    
     # attempt to find GATK4 name for POP_AF field
     # if POP_AF does not exist
    if (!POPAF.info.field %in% names(info(vcf)) &&
        "POPAF" %in% names(info(vcf))) {
        # try to add an DP geno field if missing
        POPAF.info.field <- "POPAF"
    }
    
    if (!.checkVcfFieldAvailable(vcf, "DP")) {
        # try to add an DP geno field if missing
        vcf <- .addDpField(vcf)
    }
    # check for NAs in DP
    idx <- is.na(rowSums(geno(vcf)$DP)) 
    if (any(idx)) {
        n <- .countVariants(vcf)
        vcf <- .removeVariants(vcf, idx, "NA DP")
        flog.warn("DP FORMAT field contains NAs. Removing %i variants.", n-length(vcf))
    }
    if (check.DB) {
        if (is.null(info(vcf)[[DB.info.flag]])) {
            # try to add an DB field based on rownames
            vcf <- .addDbField(vcf, DB.info.flag, POPAF.info.field, min.pop.af)
        } else if (!is.null(info(vcf)[[POPAF.info.field]])) {
            flog.warn("VCF contains both %s and %s INFO fields. Will ignore %s.",
                DB.info.flag, POPAF.info.field, POPAF.info.field)
        }
        # check for NAs in DB
        idx <- is.na(info(vcf)[[DB.info.flag]])
        if (sum(idx)) {
            flog.warn("DB INFO flag contains NAs")
            info(vcf)[[DB.info.flag]][idx] <- FALSE
        }
        cntLikelyGL <- sum(info(vcf)[[DB.info.flag]], na.rm = TRUE)
        flog.info("%i (%.1f%%) variants annotated as likely germline (%s INFO flag).",
            cntLikelyGL, cntLikelyGL/length(vcf)*100, DB.info.flag)
        if (!cntLikelyGL) {
            .stopUserError("VCF either contains no germline variants or variants are not properly annotated.")
        }
    }
    if (is.null(geno(vcf)$AD)) {
        vcf <- .addADField(vcf)
    } else {
        vcf <- .checkADField(vcf)
    }    
    if (!.checkVcfFieldAvailable(vcf, "FA")) {
        # try to add an FA geno field if missing
        vcf <- .addFaField(vcf)
    }
    idx <- !apply(geno(vcf)$FA, 1, complete.cases)
    if (any(idx)) {
        flog.warn("Found %i variants with missing allelic fraction starting with %s. Removing them.",
            sum(idx), rownames(vcf)[idx][1])
        vcf <- .removeVariants(vcf, idx, "NA FA")
    }
    meta(header(vcf))$purecnprefix <- DataFrame(
        Value = vcf.field.prefix,
        row.names = "purecnprefix"
    )
    vcf
}
.getPureCNPrefixVcf <- function(vcf) {
    prefix <- meta(header(vcf))$purecnprefix[[1]]
    prefix <- if (is.null(prefix)) "" else prefix
    return(prefix)
}
.checkADField <- function(vcf) {
    refs <- apply(geno(vcf)$AD, 2, function(x) sapply(x, function(y) y[2]))
    if (!any(complete.cases(refs))) {
        # VarScan2 does only provide alt in AD
        flog.warn("AD field misses ref counts.")
        matrixAD <- do.call(cbind, lapply(samples(header(vcf)), function(j) {
            dp <- unlist(geno(vcf)$DP[,j])
            alt <- sapply(geno(vcf)$AD[,j], function(x) x[1])
            ref <- dp - alt
            AD <- lapply(seq_along(dp), function(i) as.integer(c(ref[i], alt[i])))
            names(AD) <- names(dp)
            AD
        })) 
        colnames(matrixAD) <- samples(header(vcf))   
        geno(vcf)$AD <- matrixAD
    }
    vcf
}
    
.addDbField <- function(vcf, DB.info.flag, 
                        POPAF.info.field,
                        min.pop.af) {

    if (!is.null(info(vcf)$SOMATIC) &&
        sum(colSums(geno(vcf)$DP) > 0) == 2 &&
        any(info(vcf)$SOMATIC)) {
        db <- !info(vcf)$SOMATIC
        flog.warn("vcf.file has no DB info field for membership in germline databases.%s",
           " Found and used somatic status instead.")
    } else if (!is.null(info(vcf)[[POPAF.info.field]]) && 
        max(unlist(info(vcf)[[POPAF.info.field]]), na.rm = TRUE) > 0.05 ) {
        if (max(unlist(info(vcf)[[POPAF.info.field]]), na.rm = TRUE) > 1.1) {
            flog.info("Maximum of POPAP INFO is > 1, assuming -log10 scaled values") 
            db <- info(vcf)[[POPAF.info.field]] < -log10(min.pop.af)
        } else {    
            db <- info(vcf)[[POPAF.info.field]] > min.pop.af
        }
        db <- sapply(db, function(x) x[[1]])
        flog.warn("vcf.file has no DB info field for membership in germline databases. Found and used valid population allele frequency > %f instead.",
                  min.pop.af)
    } else {
        db <- grepl("^rs",rownames(vcf))
        
        if (!sum(db)) {
           .stopUserError("vcf.file has no %s or %s info field for membership in germline databases.", 
            DB.info.flag, POPAF.info.field)
        } else {
           flog.warn("vcf.file has no %s or %s info field for membership in germline databases.%s",
               DB.info.flag, POPAF.info.field, " Guessing it based on available dbSNP ID.")
        }
    }
    newInfo <- DataFrame(
        Number = 0, 
        Type = "Flag",
        Description = "dbSNP Membership",
        row.names = DB.info.flag)
    info(header(vcf)) <- rbind(info(header(vcf)), newInfo)
    info(vcf)[[DB.info.flag]] <- db
    vcf
}

.addSomaticField <- function(vcf) {
    # add SOMATIC flag only for GATK4 MuTect2 output
    if (ncol(vcf) < 2) {
        return(vcf)
    } else {
        newInfo <- DataFrame(
            Number=0, Type="Flag",
            Description="Somatic event",
            row.names="SOMATIC")
        info(header(vcf)) <- rbind(info(header(vcf)), newInfo)
        if (!is.null(info(vcf)$P_GERMLINE)) {
            info(vcf)$SOMATIC <- unlist(info(vcf)$P_GERMLINE < log10(0.5))
            info(vcf)$SOMATIC[is.infinite(info(vcf)$SOMATIC) | 
                is.na(info(vcf)$SOMATIC) ] <- FALSE
        } else if (!is.null(info(vcf)$GERMQ)) {
            flog.warn("Found GERMQ info field with Phred scaled germline probabilities.")
            info(vcf)$SOMATIC <- unlist(info(vcf)$GERMQ > -10* log10(0.5))
            info(vcf)$SOMATIC[is.infinite(info(vcf)$SOMATIC) | 
                is.na(info(vcf)$SOMATIC) ] <- FALSE
        } else {
            flog.warn("Having trouble guessing SOMATIC status...")
            info(vcf)$SOMATIC <- apply(geno(vcf)$GT,1,function(x) x[1]!=x[2])
        }    
    }    
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
    cosmic.vcf <- readVcf(cosmic.vcf.file, genome = genome(vcfRenamedSL)[1],  
        ScanVcfParam(which = rowRanges(vcfRenamedSL),
            info="CNT",
            fixed="ALT",
            geno=NA
        )
    )
    ov <- findOverlaps(vcfRenamedSL, cosmic.vcf, type = "equal")
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
    prefix <- .getPureCNPrefixVcf(vcf) 
    key <- paste0(prefix, "OnTarget")
    newInfo <- DataFrame(
        Number = 1, Type = "Integer",
        Description = "1: On-target; 2: Flanking region; 0: Off-target.",
        row.names = key)
    
    info(header(vcf)) <- rbind(info(header(vcf)), newInfo)
    info(vcf)[[key]] <- 0
    info(vcf)[[key]][idxPadding] <- 2
    info(vcf)[[key]][idxTarget] <- 1
    
    # report stats in log file
    targetsWithSNVs <- overlapsAny(target.granges.padding, vcf)
    percentTargetsWithSNVs <- sum(targetsWithSNVs,na.rm=TRUE)/
        length(targetsWithSNVs)*100
    flog.info("%.1f%% of targets contain variants.", 
        percentTargetsWithSNVs)
    vcf
}

.filterVcfByBQ <- function(vcf, tumor.id.in.vcf, min.base.quality) {
    n.vcf.before.filter <- .countVariants(vcf)
    idx <- NULL
    if (!is.null(geno(vcf)$BQ)) {
        # Mutect 1
        idx <- as.numeric(geno(vcf)$BQ[,tumor.id.in.vcf]) < min.base.quality
    } else if (!is.null(geno(vcf)$MBQ)) {
        # Mutect 2
        idx <- as.numeric(geno(vcf)$MBQ[,tumor.id.in.vcf]) < min.base.quality
    } else if (!is.null(rowRanges(vcf)$QUAL)) {
        # Freebayes
        idx <- as.numeric(rowRanges(vcf)$QUAL) < min.base.quality
    }
    vcf <- .removeVariants(vcf, idx, "BQ", na.rm = FALSE)
    flog.info("Removing %i low quality variants with BQ < %i.", 
        n.vcf.before.filter - .countVariants(vcf), min.base.quality) 
    vcf
}         

.checkVcfNotEmpty <- function(vcf) {
    if (!.countVariants(vcf)) .stopUserError("None of the variants in provided VCF passed filtering.")
}

# in future this will use the FILTER flag to remove variants
.removeVariants <- function(vcf, idx, label, na.rm = TRUE) {
    if (is(idx, "integer")) {
        idx <- seq(length(vcf)) %in% idx
    }    
    if (any(is.na(idx))) {
        flog.warn("Variant ids contain NAs at filter step %s.", label)
        idx[is.na(idx)] <- na.rm
    }
    # TODO: once used flags, make sure this includes already filtered variants
    if (all(idx)) {
        .stopUserError("No variants passed filter ", label, ".")
    }    
    vcf[!idx]
}
# in future this will use the FILTER flag to count
.countVariants <- function(vcf) {
    if (!is(vcf, "VCF")) {
        .stopRuntimeError("Not a VCF object in .countVariants.")
    }    
    nrow(vcf)
}    
