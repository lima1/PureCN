runAbsoluteCN <-
structure(function(# Run PureCN implementation of ABSOLUTE
### This function takes as input tumor and normal control coverage and 
### allelic fractions of germline variants and somatic mutations.
### Coverage data is provided in GATK DepthOfCoverage format, allelic fraction
### in VCF format (e.g. obtained by MuTect). Normal control does not need to 
### be matched (from the same patient). In case VCF does not contain somatic 
### status, it should contain dbSNP and optionally COSMIC annotation.
### Returns purity and ploidy combinations, sorted by likelihood score.
### Provides copy number and LOH data, by both gene and genomic region. 
gatk.normal.file=NULL, 
### GATK coverage file of normal control (optional if 
### log.ratio is provided - then it will be only used to filter low coverage 
### exons). 
### Should be already GC-normalized with \code{\link{correctCoverageBias}}. 
### Needs to be either a file name or data read with the 
### \code{\link{readCoverageGatk}} function.
##seealso<< \code{\link{correctCoverageBias} \link{segmentationCBS}}
gatk.tumor.file=NULL, 
### GATK coverage file of tumor. If NULL, requires \code{seg.file}
### and an interval file via \code{gc.gene.file}.
### Should be already GC-normalized with \code{\link{correctCoverageBias}}. 
### Needs to be either a file name or data read with the 
### \code{\link{readCoverageGatk}} function.
log.ratio=NULL, 
### Copy number log-ratios for all exons in the coverage files. 
### If NULL, calculated based on coverage files.
seg.file=NULL,
### Segmented data. Optional, to support matched SNP6 data. 
### If NULL, use coverage files or log.ratio to segment the data.  
seg.file.sdev=0.4,
### If seg.file provided, the log-ratio standard deviation, 
### used to model likelihood of sub-clonal copy number events.
vcf.file=NULL, 
### VCF file, tested with MuTect output files.  Optional, but
### typically needed to select between local optima of similar likelihood. Can 
### also be a CollapsedVCF, read with the \code{readVcf} function. 
### Requires a DB info flag for dbSNP membership. The default 
### \code{fun.setPriorVcf} function will also
### look for a Cosmic.CNT slot, containing the hits in the COSMIC database.
### Again, do not expect very useful results without a VCF file.
genome,
### Genome version, required for the \code{readVcf} function.
centromeres=NULL,
### Data.frame with centromere positions in first three columns. 
### If NULL, use pre-stored positions for genome versions hg18, hg19 and hg38. 
sex=c("?","F","M","diploid"),
### Sex of sample. If ?, detect using \code{\link{getSexFromCoverage}} 
### function and default parameters. 
### Default parameters might not work well with every assay and 
### might need to be tuned. If set to diploid, then PureCN will assume all 
### chromosomes are diploid and will not try to detect sex. 
fun.filterVcf=filterVcfMuTect, 
### Function for filtering variants. Expected output is a 
### list with elements vcf (CollapsedVCF), flag (TRUE/FALSE) and flag_comment 
### (string). The flags will be added to the output data and can be used to 
### warn users, for example when samples look too noisy. Default filter will 
### remove variants flagged by MuTect, but will keep germline variants. If 
### ran in matched normal mode, it will by default use somatic status of 
### variants and filter non-somatic calls with allelic fraction significantly 
### different from 0.5 in normal. 
args.filterVcf=list(),
### Arguments for variant filtering function. Arguments 
### vcf, tumor.id.in.vcf, coverage.cutoff and verbose are required in the 
### filter function and are automatically set (do NOT set them here again).
fun.setPriorVcf=setPriorVcf,
### Function to set prior for somatic status for each 
### variant in the VCF.
args.setPriorVcf=list(),
### Arguments for somatic prior function.
fun.segmentation=segmentationCBS, 
### Function for segmenting the copy number log-ratios. 
### Expected return value is a list with elements seg (the segmentation) and 
### size (the size in bp for all segments).
args.segmentation=list(),
### Arguments for segmentation function. Arguments 
### normal, tumor, log.ratio, plot.cnv, coverage.cutoff, sampleid, vcf, 
### tumor.id.in.vcf, verbose are required in the segmentation function and 
### automatically set (do NOT set them here again).
fun.focal=findFocal,
### Function for identifying focal amplifications.
args.focal=list(),
### Arguments for focal amplification function.
sampleid=NULL, 
### Sample id, provided in output files etc.
min.ploidy=1, 
### Minimum ploidy to be considered.
max.ploidy=6, 
### Maximum ploidy to be considered.
test.num.copy=0:7, 
### Copy numbers tested in the grid search. Note that focal 
### amplifications can have much higher copy numbers, but they will be labeled
### as subclonal (because they do not fit the integer copy numbers).
test.purity=seq(0.15,0.95,by=0.01), 
### Considered tumor purity values. 
prior.purity=rep(1,length(test.purity))/length(test.purity), 
### Priors for purity if they are available. Only change 
### when you know what you are doing.
max.candidate.solutions=20, 
### Number of local optima considered in optimization 
### and variant fitting steps. If there are too many local optima, it will use
### specified number of top candidate solutions, but will also include all 
### optima close to diploid, because silent genomes have often lots of local 
### optima.
candidates=NULL, 
### Candidates to optimize from a previous run 
### (\code{return.object$candidates}). 
### If NULL, do 2D grid search and find local optima. 
coverage.cutoff=15, 
### Minimum exon coverage in both normal and tumor. Exons
### with lower coverage are ingored. The cutoff choice depends on the expected
### purity and overall coverage. High purity samples might need a lower cutoff
### to call homozygous deletions. If an \code{exon.weight.file} (below) is NOT 
### specified, it is recommended to set a higher cutoff (e.g. 20) to remove 
### noise from unreliable exon measurements. 
max.non.clonal=0.2, 
### Maximum genomic fraction assigned to a subclonal copy 
### number state.
max.homozygous.loss=0.1,
### Maximum genomic fraction assigned to homozygous loss. 
### This is set to a fairly high default value to not exclude correct
### solutions, especially in noisy segmentations. 
iterations=30, 
### Maximum number of iterations in the Simulated Annealing copy 
### number fit optimization.
log.ratio.calibration=0.25,
### Re-calibrate log-ratios in the window 
### sd(log.ratio)*log.ratio.calibration.
remove.off.target.snvs=TRUE,
### If set to a true value, will remove all SNVs outside the 
### covered regions.
gc.gene.file=NULL, 
### A mapping file that assigns GC content and gene symbols 
### to each exon in the coverage files. Used for generating gene-level calls. 
### First column in format CHR:START-END. Second column GC content (0 to 1). 
### Third column gene symbol. This file
### can be generated with the \sQuote{GATK GCContentByInterval} tool or 
### with the \code{\link{calculateGCContentByInterval}} function.
filter.lowhigh.gc.exons=0.001,
### Quantile q (defines lower q and upper 1-q) 
### for removing exons with outlier GC profile. Assuming that GC correction 
### might not have been worked on those. Requires gc.gene.file.
max.dropout=c(0.95,1.1),
### Measures GC bias as ratio of coverage in AT-rich (GC < 0.5) 
### versus GC-rich regions (GC >= 0.5). High drop-out might indicate that
### data was not GC-normalized or that the sample quality might be 
### insufficient. Requires gc.gene.file.
filter.targeted.base=4,
### Exclude exons with targeted base (size) smaller 
### than this cutoff. This is useful when the same interval file was used to
### calculate GC content. For such small exons, the GC content is likely 
### very different from the true GC content of the probes.
max.logr.sdev=0.75,
### Flag noisy samples with segment log-ratio standard deviation 
### larger than this. Assay specific and needs to be calibrated.
max.segments=200,
### Flag noisy samples with a large number of segments. Assay 
### specific and needs to be calibrated.
prior.K=0.999,
### This defines the prior probability that the multiplicity of
### a SNV corresponds to either the maternal or the paternal copy
### number (for somatic variants additionally a multiplicity of 1).
### For perfect segmentations, this value would be 1; values smaller
### than 1 thus may provide some robustness against segmentation errors.
chr.hash=NULL,
### Mapping of non-numerical chromsome names to numerical names
### (e.g. chr1 to 1, chr2 to 2, etc.). If NULL, assume chromsomes
### are properly ordered.   
plot.cnv=TRUE, 
### Generate segmentation plots.
verbose=TRUE, 
### Verbose output.
post.optimize=FALSE,
### Optimize purity using final SCNA-fit and SNVs. This 
### might take a long time when lots of SNVs need to be fitted, but will 
### typically result in a slightly more accurate purity, especially for rather
### silent genomes or very low purities. Otherwise, it will just use the 
### purity determined via the SCNA-fit.
... 
### Additional parameters passed to the segmentation function.
) {
    debug <- FALSE

    if (missing(genome)) {
        genome <- "hg19"
        message("Default of genome=hg19 is deprecated. Please specify genome.")
    }    
    
    if (is.null(centromeres)) {
        data(centromeres, envir = environment())
        if (genome %in% names(centromeres)) {
            centromeres <- centromeres[[genome]]
        } else {
            centromeres <- NULL
        }
    }
    # argument checking
    .checkParameters(test.purity, min.ploidy, max.ploidy, max.non.clonal)
   
    test.num.copy <- sort(test.num.copy)

    if (!is.null(gc.gene.file) && !file.exists(gc.gene.file)) { 
        warning("gc.gene.file ", gc.gene.file, 
            " not found. You won't get gene-level calls.")
        gc.gene.file=NULL
    }    

    if(verbose) message("Loading GATK coverage files...")

    if (!is.null(gatk.normal.file)) {    
        if (is.character(gatk.normal.file)) {
            gatk.normal.file <- normalizePath(gatk.normal.file, mustWork=TRUE)
            normal <- readCoverageGatk(gatk.normal.file)
        } else {
            normal <- gatk.normal.file
        }    
    }
    
    if (is.null(gatk.tumor.file)) {
        if (is.null(seg.file) || is.null(gc.gene.file)) { 
            .stopUserError("Missing gatk.tumor.file requires seg.file",
                " and gc.gene.file.")
        }
        tumor <- .gcGeneToCoverage(gc.gene.file, coverage.cutoff+1)
        gatk.tumor.file <- tumor
    } else if (is.character(gatk.tumor.file)) {
        gatk.tumor.file <- normalizePath(gatk.tumor.file, mustWork=TRUE)
        tumor  <- readCoverageGatk(gatk.tumor.file)
        if (is.null(sampleid)) sampleid <- basename(gatk.tumor.file)
    } else {
        if (verbose) message(
            "gatk.tumor.file does not appear to be a filename, assuming",
            " it is valid GATK coverage data.")
        tumor <- gatk.tumor.file
        if (is.null(sampleid)) sampleid <- "Sample.1"
    }    
    if (is.null(chr.hash)) chr.hash <- .getChrHash(tumor$chr)

    # check that normal is not the same as tumor (only if no log-ratio or 
    # segmentation is provided, in that case we wouldn't use normal anyway)
    if (!is.null(gatk.normal.file) & is.null(log.ratio) & is.null(seg.file)) {
        if (identical(tumor$average.coverage, normal$average.coverage)) { 
            .stopUserError("Tumor and normal are identical. This won't give any", 
                " meaningful results and I'm stopping here.")
        }
    }    

    # this ugly if else chain covers the 3 possible ways of segmenting the 
    # data: if there is no log-ratio provided, we either need to calculate 
    # (case 1) or create fake log-ratios from a provided segmentation file 
    # (case 2). Otherwise we just take the provided log-ratio (case 3).
    if (is.null(log.ratio)) {
        if (!is.null(seg.file)) {
            if (is.null(gatk.normal.file)) normal <- tumor
            log.ratio <- .createFakeLogRatios(tumor, seg.file, chr.hash)
            if (is.null(sampleid)) sampleid <- read.delim(seg.file)[1,1]
        } else {
            if (is.null(gatk.normal.file)) {
                .stopUserError(
                "Need a normal coverage file if log.ratio and seg.file are not",
                " provided.")
            }
            log.ratio <- .calcLogRatio(normal, tumor, verbose=debug)
        }
    } else {
        # the segmentation algorithm will remove exons with low coverage in 
        # both tumor and normal, so we just use tumor if there is no normal 
        # coverage file.
        if (is.null(gatk.normal.file)) normal <- tumor
        if (!is.null(seg.file)) {
            .stopUserError("Provide either log.ratio or seg.file, not both.") 
        }
    }        
    sex <- match.arg(sex)
    sex.chr <- .getSexChr(tumor$chr)
    if (sex != "diploid") {
        if (sex =="?") {
            sex.tumor <- getSexFromCoverage(tumor, verbose=FALSE)
            sex.normal <- getSexFromCoverage(normal, verbose=FALSE)
            if (is.na(sex.tumor)) {
                tumor <- .removeChr(tumor, remove.chrs=sex.chr)
            }    
            if (is.na(sex.normal)) {
                normal <- .removeChr(normal, remove.chrs=sex.chr)
            }    
            if (!identical(sex.tumor, sex.normal)) {
                warning("Sex tumor/normal mismatch: tumor = ", sex.tumor, 
                    " normal = ", sex.normal)
            }
            sex <- sex.tumor    
            if (is.na(sex)) sex = "?"
        } 
        if (sex=="M") {
            tumor <- .removeChr(tumor, remove.chrs=sex.chr)
        } else if (sex=="F") {
            tumor <- .removeChr(tumor, remove.chrs=sex.chr[2])
        }       
    }      
          
    # NA's in log.ratio confuse the CBS function
    idx <- !is.na(log.ratio) & !is.infinite(log.ratio)
    log.ratio <- log.ratio[idx]    
    normal <- normal[idx,]
    tumor <- tumor[idx,]

    if (!is.null(gc.gene.file)) {
        gc.data <- read.delim(gc.gene.file, as.is=TRUE)
        if (!length(intersect(gc.data[,1], tumor[,1]))) {
            .stopUserError("Intervals of gatk.tumor.file and gc.gene.file ",
                "do not align.")
        }
        gc.data <- gc.data[match(as.character(tumor[,1]),
            gc.data[,1]),,drop=FALSE]
    }

    dropoutWarning <- FALSE
    # clean up noisy exons, but not if the segmentation was already provided.
    if (is.null(seg.file)) {
        if (!is.null(filter.targeted.base)) {
            idx <- which(tumor$targeted.base >= filter.targeted.base)
            if (verbose) message("Removing ", nrow(tumor)-length(idx), 
                " small exons.")
            log.ratio <- log.ratio[idx]    
            normal <- normal[idx,]
            tumor <- tumor[idx,]
        }    
        
        if (!is.null(gc.gene.file)) {
            gc.data <- gc.data[match(as.character(tumor[,1]), gc.data[,1]),]
            qq <- quantile(gc.data$gc_bias, p=c(filter.lowhigh.gc.exons, 
                1-filter.lowhigh.gc.exons), na.rm=TRUE)

            idx <- which(!(gc.data$gc_bias < qq[1] | gc.data$gc_bias > qq[2]))
            if (verbose) message("Removing ", nrow(gc.data)-length(idx),
                 " low/high GC exons.")
            gc.data <- gc.data[idx,]    
            log.ratio <- log.ratio[idx]    
            normal <- normal[idx,]
            tumor <- tumor[idx,]
            gcMetricNormal <- .calcGCmetric(gc.data, normal)
            gcMetricTumor <- .calcGCmetric(gc.data, tumor)
            if (verbose) { 
                message("AT/GC dropout: ", round(gcMetricTumor, digits=2),
                " (tumor), ", round(gcMetricNormal, digits=2), " (normal).")
            }    
            if (gcMetricNormal < max.dropout[1] || 
                gcMetricNormal > max.dropout[2] ||
                gcMetricTumor  < max.dropout[1] ||
                gcMetricTumor  > max.dropout[2]) {
                warning("High GC-bias in normal or tumor. Is data GC-normalized?")
                dropoutWarning <- TRUE
            }
        } else if (verbose) {
            message("No gc.gene.file provided. Cannot check if data was GC-normalized. Was it?")    
        }    
    }
    
    if (!is.null(gc.gene.file) && is.null(gc.data$Gene) ) {
        if (verbose) message("No Gene column in gc.gene.file.",
            " You won't get gene-level calls.")
        gc.gene.file <- NULL
    }
        
    exon.gr <- GRanges(seqnames=tumor$chr, 
        IRanges(start=tumor$probe_start,end=tumor$probe_end))

    # these are not used for segmentation, so be a little bit more careful 
    # with removing (outliers will be smoothed anyways)
    exon.well.covered.gr <- 
        exon.gr[.getWellCoveredExons(normal, tumor, coverage.cutoff*0.75)]

    vcf <- NULL
    vcf.germline <- NULL
    tumor.id.in.vcf <- NULL
    normal.id.in.vcf <- NULL
    prior.somatic <- NULL
    vcf.filtering <-list(flag=FALSE, flag_comment="")
    sex.vcf <- NULL

    if (!is.null(vcf.file)) {
        if (verbose) message("Loading VCF...")
        if (class(vcf.file) == "character") {    
            vcf <- .readAndCheckVcf(vcf.file, genome=genome)
        } else if (class(vcf.file) != "CollapsedVCF") {
            .stopUserError("vcf.file neither a filename nor a CollapsedVCF ", 
                "object.") 
        } else {
            vcf <- vcf.file
        } 

        if (length(intersect(tumor$chr,seqlevels(vcf))) < 1) {
            .stopUserError("Different chromosome names in coverage and VCF.")
        }    

        if (is.null(args.filterVcf$use.somatic.status)) {
            args.filterVcf$use.somatic.status <- TRUE
        }

        if (sum(colSums(geno(vcf)$DP)>0) == 1 && 
            args.filterVcf$use.somatic.status) { 
            message("VCF file seems to have only one sample. ", 
                "Using SNVs in single mode.")
            args.filterVcf$use.somatic.status <- FALSE
        } 

        if (args.filterVcf$use.somatic.status) { 
            normal.id.in.vcf <- names( which.max(colSums(geno(vcf)$GT=="0")) )
        }    
            
        tumor.id.in.vcf <- names( which.min(colSums(geno(vcf)$GT=="0")) )

        if (verbose) message("Assuming ", tumor.id.in.vcf, 
            " is tumor in VCF file.")
        if (sex != "diploid") {    
            sex.vcf <- getSexFromVcf(vcf, tumor.id.in.vcf, verbose=verbose)
            if (!is.na(sex.vcf) && sex %in% c("F", "M") && sex.vcf != sex) { 
                warning("Sex mismatch of coverage and VCF. ",
                "Could be because of noisy data, contamination, ",
                "loss of chrY or a mis-alignment of coverage and VCF.")
            }    
        }
        n.vcf.before.filter <- nrow(vcf)
        if (verbose) message("Found ", n.vcf.before.filter, 
            " variants in VCF file.")
        
        args.filterVcf <- c(list(vcf=vcf, tumor.id.in.vcf=tumor.id.in.vcf, 
            coverage.cutoff=coverage.cutoff, verbose=verbose), args.filterVcf)
        vcf.filtering <- do.call(fun.filterVcf, args.filterVcf)

        vcf <- vcf.filtering$vcf
        
        if (remove.off.target.snvs) {
            n.vcf.before.filter <- nrow(vcf)
            # make sure all SNVs are in covered exons
            vcf <- vcf[seq_len(nrow(vcf)) %in% queryHits(findOverlaps(vcf, 
                exon.gr))]
            if (verbose) message("Removing ", n.vcf.before.filter - nrow(vcf), 
                " variants outside intervals.", 
                " Set remove.off.target.snvs=FALSE to include.")
        }        
        args.setPriorVcf <- c(list(vcf=vcf, verbose=verbose), args.setPriorVcf) 
        prior.somatic <- do.call(fun.setPriorVcf, args.setPriorVcf)
        vcf.germline <- vcf[which(prior.somatic < 0.5)]
    }

    if (verbose) message("Sex of sample: ", sex)

    if (verbose) message("Segmenting data...")

    args.segmentation <-  c(list(normal=normal, tumor=tumor, 
        log.ratio=log.ratio, plot.cnv=plot.cnv, 
        coverage.cutoff=ifelse(is.null(seg.file), coverage.cutoff, -1), 
        sampleid=sampleid, vcf=vcf.germline, tumor.id.in.vcf=tumor.id.in.vcf,
        normal.id.in.vcf=normal.id.in.vcf, max.segments=max.segments,
        chr.hash=chr.hash, verbose=verbose,...), args.segmentation)

    vcf.germline <- NULL
    seg.result <- do.call(fun.segmentation, args.segmentation) 

    seg <- seg.result$seg
    segChrom <- seg$chrom
    seg.gr <- GRanges(seqnames=.add.chr.name(seg$chrom, chr.hash), 
        IRanges(start=round(seg$loc.start), end=seg$loc.end))
    
    snv.lr <- NULL

    if (!is.null(vcf.file)) {

        ov <- findOverlaps(seg.gr, vcf)
        n.vcf.before.filter <- nrow(vcf)
        # make sure all SNVs are in covered segments
        vcf <- vcf[subjectHits(ov)]
        args.setPriorVcf$vcf <- vcf
        args.setPriorVcf$verbose <- FALSE
        prior.somatic <- do.call(fun.setPriorVcf, args.setPriorVcf)
        if (verbose) message("Removing ", n.vcf.before.filter - nrow(vcf), 
            " variants outside segments.")
        # Add segment log-ratio to off-target snvs. 
        # For on-target, use observed log-ratio.   
        sd.ar <- sd(unlist(geno(vcf)$FA[,tumor.id.in.vcf]))
        snv.lr <- seg$seg.mean[queryHits(ov)]
        ov.vcfexon <- findOverlaps(vcf, exon.gr)
        snv.lr[queryHits(ov.vcfexon)] <- log.ratio[subjectHits(ov.vcfexon)]
        ov <- findOverlaps(seg.gr, vcf)
        if (verbose) message("Using ", nrow(vcf), " variants.")
    }
    
    # get exon log-ratios for all segments 
    ov.se <- findOverlaps(seg.gr, exon.well.covered.gr)
    exon.lrs <- lapply(seq_len(nrow(seg)), function(i) 
        log.ratio[subjectHits(ov.se)[queryHits(ov.se)==i]])
    exon.lrs <- lapply(exon.lrs, function(x) 
        subset(x, !is.na(x) & !is.infinite(x)))

    # estimate stand. dev. for exon logR within exons. this will be used as 
    # proxy for sample error.
    sd.seg <- median(sapply(exon.lrs, sd), na.rm=TRUE)
    
    # if user provided seg file, then we do not have access to the log-ratios 
    # and need to use the user provided noise estimate also, don't do outlier 
    # smoothing when we use already segmented data
    if (!is.null(seg.file)) {
        sd.seg <- seg.file.sdev
    } else {   
        exon.lrs <- lapply(exon.lrs, .smoothOutliers)
    }    

    # renormalize, in case segmentation function changed means
    exon.lrs <- lapply(seq_along(exon.lrs), function(i) {
        if (length(exon.lrs[[i]]) < 5) return(exon.lrs[[i]])
        exon.lrs[[i]]-mean(exon.lrs[[i]])+seg$seg.mean[i]
    })

    max.exon.ratio <- 7

    # show log-ratio histogram
    if (plot.cnv) { 
        if (!is.null(seg.file)) {
            seg.orig <- read.delim(seg.file)
            par(mfrow=c(2,1))
            hist(do.call(c, lapply(seq_len(nrow(seg.orig)), function(i) 
                rep(seg.orig$seg.mean[i], seg.orig$num.mark[i]))), breaks=100,
                xlab="log2 ratio", main=paste(sampleid, 
                "(original segmentation)"))
        }    
        hist(do.call(c, lapply(seq_len(nrow(seg)), function(i) 
                rep(seg$seg.mean[i], seg$num.mark[i]))), breaks=100,
                xlab="log2 ratio", main=sampleid)
        par(mfrow=c(1,1))
    }

    # initialize variables
    llik <- -Inf
    li <- seg.result$size
    C <- rep(2, length(li))

    if (sum(li < 0) > 0) .stopRuntimeError("Some segments have negative size.")

    if(verbose) { 
        message("Mean standard deviation of log-ratios: ", 
           round(sd.seg, digits=2))
    }    
    log.ratio.offset <- rep(0, nrow(seg))
     
    if(verbose) message("Optimizing purity and ploidy. ",
        "Will take a minute or two...")
    
    # find local maxima. use a coarser grid for purity, otherwise we will get
    # far too many solutions, which we will need to cluster later anyways.
    if (!is.null(candidates)){
        candidate.solutions <- candidates    
    } else {
        candidate.solutions <- .optimizeGrid(
            test.purity=seq( max(0.1,min(test.purity)), 
                min(0.9, max(test.purity)),by=1/30), min.ploidy, max.ploidy, 
            test.num.copy=test.num.copy, exon.lrs, seg, sd.seg, li, 
            max.exon.ratio, max.non.clonal, verbose, debug)

        # if we have > 20 somatic mutations, we can try estimating purity 
        # based on allelic fractions and assuming diploid genomes.
        if (!is.null(vcf.file) && sum(prior.somatic > 0.5, na.rm=TRUE) > 20) {
            somatic.purity <- min(max(test.purity), 
                .calcPuritySomaticVariants(vcf, prior.somatic, 
                    tumor.id.in.vcf))

            candidate.solutions$candidates <- rbind( 
                candidate.solutions$candidates, 
                c(2,somatic.purity ,NA,2)
            )
        }
    }    

    if (nrow(candidate.solutions$candidates) > max.candidate.solutions) {
        # test the best solutions and everything close to diploid
        idx.keep <- unique(c(seq_len(max.candidate.solutions), 
            which( ( candidate.solutions$candidates$tumor.ploidy > 1.5 &
                    candidate.solutions$candidates$tumor.ploidy < 2.6 ))))
        candidate.solutions$candidates <- 
            candidate.solutions$candidates[idx.keep,]
        warning("Too many candidate solutions! ",
            "Trying optimizing the top candidates.")
    }    
    
    if(verbose) message(paste(strwrap(paste("Local optima:", 
        paste(round(candidate.solutions$candidates$purity, digits=2),
            round(candidate.solutions$candidates$ploidy, digits=2),
            sep="/", collapse=", "))),collapse="\n"))

    simulated.annealing <- TRUE

    .optimizeSolution <- function(cpi) {

        max.attempts <- 4; attempt=0
        while (attempt < max.attempts) {
            attempt <- attempt + 1
        total.ploidy <- candidate.solutions$candidates$ploidy[cpi]
        p <- candidate.solutions$candidates$purity[cpi]

        if (verbose) message("Testing local optimum ", cpi, "/", 
            nrow(candidate.solutions$candidates),
            " at purity ", round(p, digits=2), " and total ploidy ", 
            round(total.ploidy, digits=2), ".")

        subclonal <- rep(FALSE, nrow(seg))
        old.llik <- -1; cnt.llik.equal <- 0;
        C.posterior <- matrix(ncol=length(test.num.copy)+1, nrow=nrow(seg))
        colnames(C.posterior) <- c(test.num.copy, "Subclonal")
        for (iter in seq_len(iterations)) {
            # test for convergence
            if (abs(old.llik - llik) < 0.0001) { 
                cnt.llik.equal <- cnt.llik.equal+1 
            }

            old.llik <- llik;    
            if (cnt.llik.equal > 3) break;
            subclonal.f <- length(unlist(exon.lrs[subclonal])) / 
                length(unlist(exon.lrs))
            # should not happen, but sometimes does for very unlikely local 
            # optima.
            if (subclonal.f > max.non.clonal+0.1) break
            if (iter == 1) log.ratio.offset <- .sampleOffsetFast(test.num.copy,
                 seg, exon.lrs, sd.seg, p, C, total.ploidy, max.exon.ratio, 
                 simulated.annealing, log.ratio.calibration) 

            # in the first iteration, we do not have integer copy numbers yet
            # (corresponding to local optima purity/ploidy)
            if (iter > 1) { 
                # calculate posterior probabilities of all requested purities
                total.ploidy  <- p*(sum(li*(C)))/sum(li)+(1-p)*2   #ploidy
                px.rij <- lapply(test.purity, function(px) 
                    vapply(which(!is.na(C)), function(i) 
                        .calcLlikSegment(subclonal=subclonal[i], 
                        lr=exon.lrs[[i]]+log.ratio.offset[i], sd.seg=sd.seg, 
                        p=px, Ci=C[i], total.ploidy=total.ploidy, 
                        max.exon.ratio=max.exon.ratio), double(1))
                )
                px.rij.s <- sapply(px.rij, sum, na.rm=TRUE) + log(prior.purity)
                
                if (simulated.annealing) px.rij.s <- px.rij.s * exp(iter/4)

                px.rij.s <- exp(px.rij.s-max(px.rij.s))
                # Gibbs sample purity 
                p <- test.purity[min(which(runif(n=1, min=0, 
                    max=sum(px.rij.s)) <= cumsum(px.rij.s)))]
                total.ploidy  <- p*(sum(li*(C)))/sum(li)+(1-p)*2   
                # Gibbs sample offset 
                if (iter > 2) log.ratio.offset <- .sampleOffset(subclonal, 
                    seg, exon.lrs, sd.seg, p, C, total.ploidy, max.exon.ratio, 
                    simulated.annealing, iter, log.ratio.calibration)
            }

            # calculate the log-liklihood of purity and integer copy numbers 
            # plus clonal vs subclonal status
            llik <- sum(vapply(which(!is.na(C)), function(i) 
                .calcLlikSegment(subclonal=subclonal[i], 
                    lr=exon.lrs[[i]]+log.ratio.offset[i], sd.seg=sd.seg, p=p, 
                    Ci=C[i], total.ploidy=total.ploidy, 
                    max.exon.ratio=max.exon.ratio), double(1))
            )

            if (debug) message(paste("Iteration:", iter, 
                " Log-likelihood: ", llik, 
                " Purity:",p,
                " Total Ploidy:", total.ploidy, 
                " Tumor Ploidy:", sum(li*(C))/sum(li), 
                " Fraction sub-clonal:", subclonal.f, 
                " Mean log-ratio offset", mean(log.ratio.offset)))
                 
            for (i in seq_len(nrow(seg))) { 
                # Gibbs sample copy number
                # Step 1: calculate log-likelihoods of fits
                # In the first iteration, we do not have the integer copy 
                # numbers yet, so calculate ploidy only when we have
                # it next time. Now, use the ploidy from the candidate solution. 
                if (iter > 1) total.ploidy  <- p*(sum(li*(C)))/sum(li)+(1-p)*2

                p.rij <- vapply(test.num.copy, function(Ci) 
                    .calcLlikSegment(subclonal=FALSE, 
                        lr=exon.lrs[[i]]+log.ratio.offset[i], sd.seg=sd.seg,
                        p=p, Ci=Ci, total.ploidy=total.ploidy, 
                        max.exon.ratio=max.exon.ratio), double(1))

                # calculate tumor ploidy for all possible copy numbers in this
                # segment
                ploidy <- vapply(test.num.copy, function(Ci) 
                    (sum(li[-i]*(C[-i]))+li[i]*Ci )/sum(li), double(1))

                # set probability to zero if ploidy is not within requested range 
                log.prior.ploidy <- log(ifelse(ploidy<min.ploidy | 
                    ploidy>max.ploidy,0,1))
                if (iter > 1) p.rij <- p.rij+log.prior.ploidy
                
                frac.homozygous.loss <- vapply(test.num.copy, function(Ci) 
                    (sum(li[-i] * ifelse(C[-i]==0,1,0)) + li[i] * 
                    ifelse(Ci==0,1,0))/sum(li), double(1))
                log.prior.homozygous.loss <- log(
                    ifelse(frac.homozygous.loss > max.homozygous.loss,0,1))
                if (iter > 1) p.rij <- p.rij+log.prior.homozygous.loss

                # model sub clonal state with a uniform distribution
                p.rij <- c(p.rij, .calcLlikSegmentSubClonal( 
                    exon.lrs[[i]]+log.ratio.offset[i], max.exon.ratio))

                C.posterior[i,] <- exp(p.rij-max(p.rij))  

                if (simulated.annealing) p.rij <- p.rij * exp(iter/4)
                
                # Because we are in log-space, sample relative to most likely fit 
                p.rij <- exp(p.rij-max(p.rij))
                # Now Gibbs sample best fit
                z <- runif(n=1, min=0, max=sum(p.rij))
                if (is.na(z)) {
                    message(paste(iter, i, ploidy))
                }    
                id <- min(which(z <= cumsum(p.rij)))
                old.C <- C[i]
                opt.C <- max((2^(seg$seg.mean[i])*total.ploidy)/ 
                    p-((2*(1-p))/p),0)
                if (id > length(test.num.copy)) {
                    # optimal non-integer copy number
                    C[i] <- opt.C
                    subclonal[i] <- TRUE
                } else {    
                    C[i] <- test.num.copy[id]
                    subclonal[i] <- FALSE
                }
                if (old.C != C[i] && debug) message("Old: ", old.C, " New: ", 
                    C[i], " LR: ", mean(exon.lrs[[i]]))
            }
        }
         if (subclonal.f < max.non.clonal && abs(total.ploidy - 
            candidate.solutions$candidates$ploidy[cpi]) < 1) break
         log.ratio.calibration <- log.ratio.calibration + 0.25    
         if (verbose && attempt < max.attempts) {
            message("Recalibrating log-ratios...")
         }   
        }    
        seg.adjusted <- seg
        seg.adjusted$C <- C
        seg.adjusted$size <- li
        llik.snv <- NULL
        SNV.posterior <- NULL

        if (subclonal.f > max.non.clonal) {
            if (!is.null(vcf.file)) {
                # we skipped the SNV fitting for this rare corner case.
                SNV.posterior <- list(
                    beta.model=list(llik=-Inf))
            }    
            return(list(
                log.likelihood=llik, 
                purity=p, 
                ploidy=weighted.mean(C,li), 
                total.ploidy=total.ploidy, 
                seg=seg.adjusted, 
                C.posterior=data.frame(C.posterior/rowSums(C.posterior), 
                ML.C=C, 
                ML.Subclonal=subclonal), 
                SNV.posterior=SNV.posterior, 
                fraction.subclonal=subclonal.f, 
                gene.calls=NA, 
                log.ratio.offset=log.ratio.offset))
        }

        if (!is.null(gc.gene.file)) {
            gene.calls <- .getGeneCalls(seg.adjusted, gc.data, log.ratio, 
                fun.focal, args.focal, chr.hash)
        } else {
            gene.calls <- NA
        }       

        if (!is.null(vcf.file)) {
            if (post.optimize) { 
                idx <- (max(1, 
                    match(p, test.purity)-4)):(min(length(test.purity),  
                    match(p, test.purity)+4))
                tp <- test.purity[idx]
                pp <- prior.purity[idx]
            } else {
                tp <- p
                pp <- 1
            }
            res.snvllik <- lapply(tp, function(px) {
                if (verbose) message("Fitting SNVs for purity ", 
                    round(px, digits=2), " and tumor ploidy ", 
                    round( weighted.mean(C,li),digits=2), ".")

                list(
                    beta.model  = .calcSNVLLik(vcf, tumor.id.in.vcf, ov, px, 
                        test.num.copy, C.posterior, C, snv.model="beta", 
                        prior.somatic, snv.lr, sampleid, prior.K=prior.K)
                )})

            if (post.optimize) {
                px.rij <- lapply(tp, function(px) vapply(which(!is.na(C)), 
                    function(i) .calcLlikSegment(subclonal=subclonal[i], 
                    lr=exon.lrs[[i]]+log.ratio.offset[i], sd.seg=sd.seg, p=px, 
                    Ci=C[i], total.ploidy=px*(sum(li*(C)))/sum(li)+(1-px)*2, 
                    max.exon.ratio=max.exon.ratio), double(1))
                )

                px.rij.s <- sapply(px.rij, sum, na.rm=TRUE) + log(pp) + 
                    vapply(res.snvllik, function(x) x$beta.model$llik, 
                        double(1))
            
                px.rij.s <- exp(px.rij.s-max(px.rij.s))
                idx <- which.max(px.rij.s)
            } else {
                idx <- 1
            }    
            p <- tp[idx]
            if (verbose) message("Optimized purity: ", p)
            SNV.posterior <- res.snvllik[[idx]]
        }

        list(
            log.likelihood=llik, 
            purity=p, 
            ploidy=weighted.mean(C,li), 
            total.ploidy=total.ploidy, 
            seg=seg.adjusted, 
            C.posterior=data.frame(C.posterior/rowSums(C.posterior), 
            ML.C=C, 
            ML.Subclonal=subclonal), 
            SNV.posterior=SNV.posterior, 
            fraction.subclonal=subclonal.f, 
            fraction.homozygous.loss=sum(li[which(C<0.01)])/sum(li), 
            gene.calls=gene.calls, 
            log.ratio.offset=log.ratio.offset, 
            failed=FALSE)
    }

    results <- lapply(seq_len(nrow(candidate.solutions$candidates)),
        .optimizeSolution)
    if (verbose) message(
        "Remember, posterior probabilities assume a correct SCNA fit.")

    results <- .rankResults(results)
    results <- .filterDuplicatedResults(results)
    results <- .flagResults(results, max.non.clonal=max.non.clonal, 
        max.logr.sdev=max.logr.sdev, logr.sdev=sd.seg, 
        max.segments=max.segments, flag=vcf.filtering$flag, 
        flag_comment=vcf.filtering$flag_comment, dropout=dropoutWarning)  

    if (!is.null(gc.gene.file)) {
        # Add LOH calls to gene level calls 
        for (i in seq_along(results)) {
            results[[i]]$gene.calls <- .getGeneCallsLOH( results[[i]] )
        }
    }    

    if (!is.null(vcf.file)) {
        for (i in seq_along(results)) {
            results[[i]]$SNV.posterior$beta.model$posteriors <- 
                .getVariantCallsLOH( results[[i]] )
        }
    }
    ##value<< A list with elements
    list(
        candidates=candidate.solutions, ##<< Results of the grid search.
        results=results, ##<< All local optima, sorted by final rank.
        input=list( ##<< The input data.
            tumor=gatk.tumor.file, 
            normal=gatk.normal.file, 
            log.ratio=data.frame(probe=normal[,1], log.ratio=log.ratio), 
            log.ratio.sdev=sd.seg, 
            vcf=vcf, 
            sampleid=sampleid, 
            sex=sex, 
            sex.vcf=sex.vcf,
            chr.hash=chr.hash,
            centromeres=centromeres) 
        )
##end<<
},ex=function(){
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")
vcf.file <- system.file("extdata", "example_vcf.vcf", 
    package="PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
    package="PureCN")

# Speed-up the runAbsoluteCN call by using the stored grid-search 
# (purecn.example.output$candidates).
data(purecn.example.output)

# The max.candidate.solutions parameter is set to a very low value only to
# speed-up this example.  This is not a good idea for real samples.
# We also remove off target snvs, also only to speed-up the example.

ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
    gatk.tumor.file=gatk.tumor.file, remove.off.target.snvs=TRUE,
    genome="hg19", candidates=purecn.example.output$candidates, 
    max.candidate.solutions=2, vcf.file=vcf.file, sampleid='Sample1', 
    gc.gene.file=gc.gene.file)
})    
