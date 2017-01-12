#' Run PureCN implementation of ABSOLUTE
#' 
#' This function takes as input tumor and normal control coverage and allelic
#' fractions of germline variants and somatic mutations. Coverage data is
#' provided in GATK DepthOfCoverage format, allelic fraction in VCF format
#' (e.g. obtained by MuTect). Normal control does not need to be matched (from
#' the same patient). In case VCF does not contain somatic status, it should
#' contain dbSNP and optionally COSMIC annotation. Returns purity and ploidy
#' combinations, sorted by likelihood score. Provides copy number and LOH data,
#' by both gene and genomic region.
#' 
#' 
#' @param normal.coverage.file GATK coverage file of normal control (optional
#' if log.ratio is provided - then it will be only used to filter low coverage
#' exons).  Should be already GC-normalized with
#' \code{\link{correctCoverageBias}}.  Needs to be either a file name or data
#' read with the \code{\link{readCoverageGatk}} function.
#' @param tumor.coverage.file GATK coverage file of tumor. If \code{NULL},
#' requires \code{seg.file} and an interval file via \code{gc.gene.file}.
#' Should be already GC-normalized with \code{\link{correctCoverageBias}}.
#' Needs to be either a file name or data read with the
#' \code{\link{readCoverageGatk}} function.
#' @param log.ratio Copy number log-ratios for all exons in the coverage files.
#' If \code{NULL}, calculated based on coverage files.
#' @param seg.file Segmented data. Optional, to support matched SNP6 data.  If
#' \code{NULL}, use coverage files or \code{log.ratio} to segment the data.
#' @param seg.file.sdev If \code{seg.file} provided, the log-ratio standard
#' deviation, used to model likelihood of sub-clonal copy number events.
#' @param vcf.file VCF file, tested with \sQuote{MuTect 1} output files.
#' Optional, but typically needed to select between local optima of similar
#' likelihood. Can also be a \code{CollapsedVCF}, read with the \code{readVcf}
#' function.  Requires a DB info flag for dbSNP membership. The default
#' \code{fun.setPriorVcf} function will also look for a Cosmic.CNT slot (see
#' \code{cosmic.vcf.file}), containing the hits in the COSMIC database. Again,
#' do not expect very useful results without a VCF file.
#' @param normalDB Normal database, created with
#' \code{\link{createNormalDatabase}}. If provided, used to calculate gene-level
#' p-values (requires \code{Gene} column in \code{gc.gene.file}) and to filter
#' targets with low coverage in the pool of normal samples.
#' @param genome Genome version, for example hg19.
#' @param centromeres A \code{data.frame} with centromere positions in first
#' three columns.  If \code{NULL}, use pre-stored positions for genome versions
#' hg18, hg19 and hg38.
#' @param sex Sex of sample. If \code{?}, detect using
#' \code{\link{getSexFromCoverage}} function and default parameters.  Default
#' parameters might not work well with every assay and might need to be tuned.
#' If set to diploid, then PureCN will assume all chromosomes are diploid and
#' will not try to detect sex.
#' @param fun.filterVcf Function for filtering variants. Expected output is a
#' list with elements \code{vcf} (\code{CollapsedVCF}), flag
#' (\code{logical(1)}) and \code{flag_comment} (\code{character(1)}). The flags
#' will be added to the output data and can be used to warn users, for example
#' when samples look too noisy. Default filter will remove variants flagged by
#' MuTect, but will keep germline variants. If ran in matched normal mode, it
#' will by default use somatic status of variants and filter non-somatic calls
#' with allelic fraction significantly different from 0.5 in normal. Defaults
#' to \code{\link{filterVcfMuTect}}, which in turn also calls
#' \code{\link{filterVcfBasic}}.
#' @param args.filterVcf Arguments for variant filtering function. Arguments
#' \code{vcf}, \code{tumor.id.in.vcf}, \code{min.coverage},
#' \code{model.homozygous} and \code{error} are required in the
#' filter function and are automatically set.
#' @param fun.setPriorVcf Function to set prior for somatic status for each
#' variant in the VCF. Defaults to \code{\link{setPriorVcf}}.
#' @param args.setPriorVcf Arguments for somatic prior function.
#' @param fun.setMappingBiasVcf Function to set mapping bias for each variant
#' in the VCF. Defaults to \code{\link{setMappingBiasVcf}}.
#' @param args.setMappingBiasVcf Arguments for mapping bias function.
#' @param fun.filterTargets Function for filtering low-quality targets in the
#' coverage files. Needs to return a \code{logical} vector whether an interval
#' should be used for segmentation. Defaults to \code{\link{filterTargets}}.
#' @param args.filterTargets Arguments for target filtering function. Arguments
#' \code{log.ratio}, \code{tumor}, \code{gc.data}, \code{seg.file} and
#' \code{normalDB} are required and automatically set 
#' @param fun.segmentation Function for segmenting the copy number log-ratios.
#' Expected return value is a \code{data.frame} representation of the
#' segmentation. Defaults to \code{\link{segmentationCBS}}.
#' @param args.segmentation Arguments for segmentation function. Arguments
#' \code{normal}, \code{tumor}, \code{log.ratio}, \code{plot.cnv} and
#' \code{min.coverage}, \code{sampleid}, \code{vcf}, \code{tumor.id.in.vcf},
#' \code{centromeres} are required in the segmentation function 
#' and automatically set.
#' @param fun.focal Function for identifying focal amplifications. Defaults to
#' \code{\link{findFocal}}.
#' @param args.focal Arguments for focal amplification function.
#' @param sampleid Sample id, provided in output files etc.
#' @param min.ploidy Minimum ploidy to be considered.
#' @param max.ploidy Maximum ploidy to be considered.
#' @param test.num.copy Copy numbers tested in the grid search. Note that focal
#' amplifications can have much higher copy numbers, but they will be labeled
#' as subclonal (because they do not fit the integer copy numbers).
#' @param test.purity Considered tumor purity values.
#' @param prior.purity \code{numeric(length(test.purity))} with priors for
#' tested purity values. If \code{NULL}, use flat priors.
#' @param prior.K This defines the prior probability that the multiplicity of a
#' SNV corresponds to either the maternal or the paternal copy number (for
#' somatic variants additionally to a multiplicity of 1). For perfect
#' segmentations, this value would be 1; values smaller than 1 thus may provide
#' some robustness against segmentation errors.
#' @param prior.contamination The prior probability that a known SNP is from a
#' different individual.
#' @param max.candidate.solutions Number of local optima considered in
#' optimization and variant fitting steps. If there are too many local optima,
#' it will use specified number of top candidate solutions, but will also
#' include all optima close to diploid, because silent genomes have often lots
#' of local optima.
#' @param candidates Candidates to optimize from a previous run
#' (\code{return.object$candidates}).  If \code{NULL}, do 2D grid search and
#' find local optima.
#' @param min.coverage Minimum coverage in both normal and tumor. Targets with
#' lower coverage are ignored.
#' @param max.coverage.vcf This will set the maximum number of reads in the SNV
#' fitting.  This is to avoid that small non-reference biases that come
#' apparent only at high coverages have a dramatic influence on likelihood
#' scores.
#' @param max.non.clonal Maximum genomic fraction assigned to a subclonal copy
#' number state.
#' @param max.homozygous.loss \code{double(2)} with maximum genomic fraction 
#' assigned to homozygous loss and maximum size of a homozygous loss segment.  
#' These are set to a fairly high default value to not exclude correct
#' solutions, especially in noisy segmentations.
#' @param non.clonal.M Average expected cellular fraction of sub-clonal somatic
#' mutations. This is to calculate expected allelic fractions of a single
#' sub-clonal bin for SNVs. For all somatic variants, more accurate cellular
#' fractions are calculated.
#' @param max.mapping.bias Exclude variants with high mapping bias from the
#' likelihood score calculation. Note that bias is reported on an inverse
#' scale; a variant with mapping bias of 1 has no bias.
#' @param iterations Maximum number of iterations in the Simulated Annealing
#' copy number fit optimization. Note that this an integer optimization problem
#' that should converge quickly. Allowed range is 10 to 250.
#' @param log.ratio.calibration Re-calibrate log-ratios in the window
#' \code{sd(log.ratio)*log.ratio.calibration}.
#' @param remove.off.target.snvs Deprecated. Use the corresponding argument in
#' \code{args.filterVcf}.
#' @param smooth.log.ratio Smooth \code{log.ratio} using the \code{DNAcopy}
#' package.
#' @param model.homozygous Homozygous germline SNPs are uninformative and by
#' default removed. In 100 percent pure samples such as cell lines, however,
#' heterozygous germline SNPs appear homozygous in case of LOH. Setting this
#' parameter to \code{TRUE} will keep homozygous SNPs and include a homozygous
#' SNP state in the likelihood model. Not necessary when matched normal samples
#' are available.
#' @param error Estimated sequencing error rate. Used to calculate minimum
#' number of supporting reads for SNVs using
#' \code{\link{calculatePowerDetectSomatic}}. Also used to calculate the
#' probability of homozygous SNP allelic fractions (assuming reference reads
#' are sequencing errors).
#' @param gc.gene.file A mapping file that assigns GC content and gene symbols
#' to each exon in the coverage files. Used for generating gene-level calls.
#' First column in format CHR:START-END. Second column GC content (0 to 1).
#' Third column gene symbol. This file can be generated with the \sQuote{GATK
#' GCContentByInterval} tool or with the
#' \code{\link{calculateGCContentByInterval}} function.
#' @param max.dropout Measures GC bias as ratio of coverage in AT-rich (GC <
#' 0.5) versus GC-rich regions (GC >= 0.5). High drop-out might indicate that
#' data was not GC-normalized or that the sample quality might be insufficient.
#' Requires \code{gc.gene.file}.
#' @param max.logr.sdev Flag noisy samples with segment log-ratio standard
#' deviation larger than this. Assay specific and needs to be calibrated.
#' @param max.segments Flag noisy samples with a large number of segments.
#' Assay specific and needs to be calibrated.
#' @param min.gof Flag purity/ploidy solutions with poor fit.
#' @param plot.cnv Generate segmentation plots.
#' @param cosmic.vcf.file Add a \code{Cosmic.CNT} info field to the provided
#' \code{vcf.file} using a VCF file containing the COSMIC database. The default
#' \code{fun.setPriorVcf} function will give SNVs found in the COSMIC database
#' a higher prior probability of being somatic. Not used in likelhood model
#' when matched normal is available in \code{vcf.file}. Should be compressed
#' and indexed with bgzip and tabix, respectively.
#' @param post.optimize Optimize purity using final SCNA-fit and SNVs. This
#' might take a long time when lots of SNVs need to be fitted, but will
#' typically result in a slightly more accurate purity, especially for rather
#' silent genomes or very low purities. Otherwise, it will just use the purity
#' determined via the SCNA-fit.
#' @param log.file If not \code{NULL}, store verbose output to file.
#' @param verbose Verbose output.
#' @return A list with elements \item{candidates}{Results of the grid search.}
#' \item{results}{All local optima, sorted by final rank.} \item{input}{The
#' input data.}
#' @author Markus Riester
#' @references Riester et al. (2016). PureCN: Copy number calling and SNV 
#' classification using targeted short read sequencing. Source Code for Biology 
#' and Medicine, 11, pp. 13. 
#'
#' Carter et al. (2012), Absolute quantification of somatic DNA alterations in 
#' human cancer. Nature Biotechnology.
#'
#' @seealso \code{\link{correctCoverageBias}} \code{\link{segmentationCBS}}
#' \code{\link{calculatePowerDetectSomatic}}
#' @examples
#' 
#' normal.coverage.file <- system.file('extdata', 'example_normal.txt', 
#'     package='PureCN')
#' tumor.coverage.file <- system.file('extdata', 'example_tumor.txt', 
#'     package='PureCN')
#' vcf.file <- system.file('extdata', 'example_vcf.vcf', 
#'     package='PureCN')
#' gc.gene.file <- system.file('extdata', 'example_gc.gene.file.txt', 
#'     package='PureCN')
#' 
#' # The max.candidate.solutions, max.ploidy and test.purity parameters are set to
#' # non-default values to speed-up this example.  This is not a good idea for real
#' # samples.
#' ret <-runAbsoluteCN(normal.coverage.file=normal.coverage.file, 
#'     tumor.coverage.file=tumor.coverage.file, genome='hg19', vcf.file=vcf.file,
#'     sampleid='Sample1', gc.gene.file=gc.gene.file,
#'     max.ploidy=4, test.purity=seq(0.3,0.7,by=0.05), max.candidate.solutions=1)
#' 
#' 
#' # If a high-quality segmentation was obtained with third-party tools:
#' seg.file <- system.file('extdata', 'example_seg.txt', 
#'     package = 'PureCN')
#' 
#' # By default, PureCN will re-segment the data, for example to identify
#' # regions of copy number neutral LOH. If this is not wanted, we can provide
#' # a minimal segmentation function which just returns the provided one:
#' funSeg <- function(seg, ...) return(seg)
#' 
#' res <- runAbsoluteCN(seg.file=seg.file, fun.segmentation=funSeg, max.ploidy = 4,
#'     test.purity = seq(0.3, 0.7, by = 0.05), max.candidate.solutions=1,
#'     genome='hg19', gc.gene.file=gc.gene.file)
#' 
#' @export runAbsoluteCN
#' @import DNAcopy
#' @import IRanges
#' @import VariantAnnotation
#' @importFrom GenomicRanges GRanges
#' @importFrom stats complete.cases dbeta dnorm dunif runif weighted.mean
#'             dbinom C
#' @importFrom utils data read.delim tail packageVersion
#' @importFrom S4Vectors queryHits subjectHits DataFrame
#' @importFrom data.table data.table
#' @importFrom futile.logger flog.info flog.warn flog.fatal
#'             flog.threshold flog.appender appender.tee
runAbsoluteCN <- function(normal.coverage.file = NULL, 
    tumor.coverage.file = NULL, log.ratio = NULL, seg.file = NULL, 
    seg.file.sdev = 0.4, vcf.file = NULL, normalDB = NULL, genome, 
    centromeres = NULL, sex = c("?", "F", "M", "diploid"), 
    fun.filterVcf = filterVcfMuTect, args.filterVcf = list(), 
    fun.setPriorVcf = setPriorVcf, args.setPriorVcf = list(), 
    fun.setMappingBiasVcf = setMappingBiasVcf, args.setMappingBiasVcf = list(), 
    fun.filterTargets = filterTargets, args.filterTargets = list(), 
    fun.segmentation = segmentationCBS, args.segmentation = list(), 
    fun.focal = findFocal, args.focal = list(), 
    sampleid = NULL, min.ploidy = 1, max.ploidy = 6, test.num.copy = 0:7, 
    test.purity = seq(0.15, 0.95, by = 0.01), prior.purity = NULL, 
    prior.K = 0.999, prior.contamination = 0.01, max.candidate.solutions = 20,
    candidates = NULL, min.coverage = 15, max.coverage.vcf = 300, 
    max.non.clonal = 0.2, max.homozygous.loss = c(0.05, 1e07) , non.clonal.M = 1/3, 
    max.mapping.bias = 0.8, iterations = 30, log.ratio.calibration = 0.25, 
    smooth.log.ratio = TRUE, remove.off.target.snvs = NULL, 
    model.homozygous = FALSE, error = 0.001, gc.gene.file = NULL, 
    max.dropout = c(0.95, 1.1), max.logr.sdev = 0.75, 
    max.segments = 300, min.gof = 0.8, plot.cnv = TRUE, cosmic.vcf.file = NULL, 
    post.optimize = FALSE, log.file = NULL, verbose = TRUE) {

    if (!verbose) flog.threshold("WARN")
    if (!is.null(log.file)) flog.appender(appender.tee(log.file))

     # log function arguments     
    try(.logHeader(as.list(match.call())[-1]), silent=TRUE)

    # TODO: remove in PureCN 1.8
    if (!is.null(remove.off.target.snvs)) {
        args.filterVcf$remove.off.target.snvs <- remove.off.target.snvs
        flog.warn("remove.off.target.snvs is deprecated. Please use it in args.filterVcf instead.")
    }
    # TODO: remove in PureCN 1.8
    if (length(max.homozygous.loss)==1){
         max.homozygous.loss <- c(max.homozygous.loss, 1e07)
         flog.warn("max.homozygous.loss now a double(2) vector. Please provide both values.")
    }     
    # TODO: remove in PureCN 1.8
    if (is.null(normalDB) && !is.null(args.filterTargets$normalDB)) {
        normalDB <- args.filterTargets$normalDB
        flog.warn("normalDB now a runAbsoluteCN argument. Please provide it there, not in args.filterTargets.")
    }    
    centromeres <- .getCentromerePositions(centromeres, genome)
    
    # defaults to equal priors for all tested purity values
    if (is.null(prior.purity)) {
        prior.purity <- rep(1, length(test.purity))/length(test.purity)
    }
    # argument checking
    .checkParameters(test.purity, min.ploidy, max.ploidy, max.non.clonal,
        max.homozygous.loss, sampleid, prior.K, prior.contamination, prior.purity,
        iterations, min.gof, model.homozygous, gc.gene.file)
    
    test.num.copy <- sort(test.num.copy)
    
    flog.info("Loading GATK coverage files...")
    
    if (!is.null(normal.coverage.file)) {
        if (is.character(normal.coverage.file)) {
            normal.coverage.file <- normalizePath(normal.coverage.file, mustWork = TRUE)
            normal <- readCoverageGatk(normal.coverage.file)
        } else {
            normal <- normal.coverage.file
        }
    }
    if (is.null(tumor.coverage.file)) {
        if ((is.null(seg.file) && is.null(log.ratio)) || is.null(gc.gene.file)) {
            .stopUserError("Missing tumor.coverage.file requires seg.file or ", 
                           "log.ratio and gc.gene.file.")
        }
        tumor <- .gcGeneToCoverage(gc.gene.file, min.coverage + 1)
        tumor.coverage.file <- tumor
    } else if (is.character(tumor.coverage.file)) {
        tumor.coverage.file <- normalizePath(tumor.coverage.file, mustWork = TRUE)
        tumor <- readCoverageGatk(tumor.coverage.file)
        if (is.null(sampleid)) 
            sampleid <- basename(tumor.coverage.file)
    } else {
        tumor <- tumor.coverage.file
        if (is.null(sampleid)) 
            sampleid <- "Sample.1"
    }
    
    chr.hash <- .getChrHash(tumor$chr)
    
    # check that normal is not the same as tumor (only if no log-ratio or
    # segmentation is provided, in that case we wouldn't use normal anyway)
    if (!is.null(normal.coverage.file) & is.null(log.ratio) & is.null(seg.file)) {
        if (identical(tumor$average.coverage, normal$average.coverage)) {
            .stopUserError("Tumor and normal are identical. This won't give any", 
                " meaningful results and I'm stopping here.")
        }
    }
    
    # this ugly if else chain covers the 3 possible ways of segmenting the data: if
    # there is no log-ratio provided, we either need to calculate (case 1) or create
    # fake log-ratios from a provided segmentation file (case 2). Otherwise we just
    # take the provided log-ratio (case 3).
    if (is.null(log.ratio)) {
        if (!is.null(seg.file)) {
            if (is.null(normal.coverage.file)) 
                normal <- tumor
            log.ratio <- .createFakeLogRatios(tumor, seg.file, chr.hash)
            smooth.log.ratio <- FALSE
            if (is.null(sampleid)) 
                sampleid <- read.delim(seg.file, as.is=TRUE)[1, 1]
        } else {
            if (is.null(normal.coverage.file)) {
                .stopUserError("Need a normal coverage file if log.ratio and seg.file are not", 
                  " provided.")
            }
            log.ratio <- calculateLogRatio(normal, tumor)
        }
    } else {
        # the segmentation algorithm will remove targets with low coverage in both tumor
        # and normal, so we just use tumor if there is no normal coverage file.
        if (is.null(normal.coverage.file)) 
            normal <- tumor
        if (length(log.ratio) != nrow(tumor)) {
            .stopUserError("Length of log.ratio different from tumor ", 
                           "coverage.")
        }
        if (is.null(sampleid)) 
            sampleid <- "Sample.1"
    }
    
    sex <- .getSex(match.arg(sex), normal, tumor)
    tumor <- .fixAllosomeCoverage(sex, tumor)
    
    gc.data <- NULL
    if (!is.null(gc.gene.file)) {
        gc.data <- read.delim(gc.gene.file, as.is = TRUE)
        if (!length(intersect(gc.data[, 1], tumor[, 1]))) {
            .stopUserError("Intervals of tumor.coverage.file and gc.gene.file ", 
                "do not align.")
        }
        gc.data <- gc.data[match(as.character(tumor[, 1]), gc.data[, 1]), , drop = FALSE]
    }
    
    args.filterTargets <- c(list(log.ratio = log.ratio, tumor = tumor, 
        gc.data = gc.data, seg.file = seg.file, normalDB = normalDB),
        args.filterTargets)
    
    targetsUsed <- do.call(fun.filterTargets, 
        .checkArgs(args.filterTargets, "filterTargets"))
    
    # chr.hash is an internal data structure, so we need to do this separately.
    targetsUsed <- .filterTargetsChrHash(targetsUsed, tumor, chr.hash)
    targetsUsed <- which(targetsUsed)
    if (nrow(tumor) != nrow(normal) || nrow(tumor) != length(log.ratio) || (!is.null(gc.gene.file) && 
        nrow(tumor) != nrow(gc.data))) {
        .stopRuntimeError("Mis-aligned coverage data.")
    }
    tumor <- tumor[targetsUsed, ]
    normal <- normal[targetsUsed, ]
    log.ratio <- log.ratio[targetsUsed]
    if (!is.null(gc.gene.file)) {
        gc.data <- gc.data[targetsUsed, ]
    }
    flog.info("Using %i targets.", nrow(tumor))

    if (smooth.log.ratio) {
        CNA.obj <- smooth.CNA(CNA(log.ratio, 
            .strip.chr.name(normal$chr, chr.hash), 
            floor((normal$probe_start + normal$probe_end)/2), 
            data.type="logratio", sampleid="sample"))
        log.ratio <- CNA.obj$sample
    }     
    
    dropoutWarning <- FALSE
    # clean up noisy targets, but not if the segmentation was already provided.
    if (is.null(seg.file)) {
        if (!is.null(gc.gene.file)) {
            dropoutWarning <- .checkGCBias(normal, tumor, gc.data, max.dropout)
        } else {
            flog.info("No gc.gene.file provided. Cannot check if data was GC-normalized. Was it?")
        }
    }
    
    if (!is.null(gc.gene.file) && is.null(gc.data$Gene)) {
        flog.info("No Gene column in gc.gene.file. You won't get gene-level calls.")
        gc.gene.file <- NULL
    }
    
    exon.gr <- GRanges(seqnames=tumor$chr, IRanges(start=tumor$probe_start, 
                                                   end=tumor$probe_end))
    
    vcf <- NULL
    vcf.germline <- NULL
    tumor.id.in.vcf <- NULL
    normal.id.in.vcf <- NULL
    prior.somatic <- NULL
    mapping.bias <- NULL
    vcf.filtering <- list(flag = FALSE, flag_comment = "")
    sex.vcf <- NULL
    
    if (!is.null(vcf.file)) {
        flog.info("Loading VCF...")
        vcf <- .readAndCheckVcf(vcf.file, genome=genome)
        
        if (length(intersect(tumor$chr, seqlevels(vcf))) < 1) {
            .stopUserError("Different chromosome names in coverage and VCF.")
        }
        
        if (is.null(args.filterVcf$use.somatic.status)) {
            args.filterVcf$use.somatic.status <- TRUE
        }
        
        if (sum(colSums(geno(vcf)$DP) > 0) == 1 && args.filterVcf$use.somatic.status) {
            args.filterVcf$use.somatic.status <- FALSE
        }
        
        tumor.id.in.vcf <- .getTumorIdInVcf(vcf, sampleid = sampleid)
        
        if (args.filterVcf$use.somatic.status) {
            normal.id.in.vcf <- .getNormalIdInVcf(vcf, tumor.id.in.vcf)
        }
        
        flog.info("%s is tumor in VCF file.", tumor.id.in.vcf)
        if (sex != "diploid") {
            sex.vcf <- getSexFromVcf(vcf, tumor.id.in.vcf)
            if (!is.na(sex.vcf) && sex %in% c("F", "M") && sex.vcf != sex) {
                flog.warn("Sex mismatch of coverage and VCF. %s%s", 
                        "Could be because of noisy data, contamination, ", 
                  "loss of chrY or a mis-alignment of coverage and VCF.")
            }
        }
        n.vcf.before.filter <- nrow(vcf)
        flog.info("Found %i variants in VCF file.", n.vcf.before.filter)
        
        args.filterVcf <- c(list(vcf = vcf, tumor.id.in.vcf = tumor.id.in.vcf, 
            model.homozygous = model.homozygous, error = error, 
            target.granges = exon.gr), args.filterVcf)
        if (is.null(args.filterVcf$min.coverage)) {
            args.filterVcf$min.coverage <- min.coverage
        }
        vcf.filtering <- do.call(fun.filterVcf, 
            .checkArgs(args.filterVcf, "filterVcf"))
        
        vcf <- vcf.filtering$vcf
        
        if (!is.null(cosmic.vcf.file)) {
            vcf <- .addCosmicCNT(vcf, cosmic.vcf.file)
        }
        
        args.setPriorVcf <- c(list(vcf = vcf, tumor.id.in.vcf = tumor.id.in.vcf),
            args.setPriorVcf)
        prior.somatic <- do.call(fun.setPriorVcf, 
            .checkArgs(args.setPriorVcf, "setPriorVcf"))
        
        # get mapping bias
        args.setMappingBiasVcf$vcf <- vcf
        args.setMappingBiasVcf$tumor.id.in.vcf <- tumor.id.in.vcf
        mapping.bias <- do.call(fun.setMappingBiasVcf,
            .checkArgs(args.setMappingBiasVcf, "setMappingBiasVcf"))
        idxHqGermline <- prior.somatic < 0.5 & mapping.bias >= max.mapping.bias
        vcf.germline <- vcf[idxHqGermline]
    }
    
    flog.info("Sample sex: %s", sex)
    flog.info("Segmenting data...")
    
    args.segmentation <- c(list(normal = normal, tumor = tumor, log.ratio = log.ratio, 
        seg = .loadSegFile(seg.file), plot.cnv = plot.cnv, min.coverage = ifelse(is.null(seg.file), 
            min.coverage, -1), sampleid = sampleid, vcf = vcf.germline, tumor.id.in.vcf = tumor.id.in.vcf, 
        normal.id.in.vcf = normal.id.in.vcf, max.segments = max.segments, chr.hash = chr.hash, 
        centromeres = centromeres), args.segmentation)
    
    vcf.germline <- NULL
    seg <- do.call(fun.segmentation,
            .checkArgs(args.segmentation, "segmentation"))
    
    seg.gr <- GRanges(seqnames = .add.chr.name(seg$chrom, chr.hash), IRanges(start = round(seg$loc.start), 
        end = seg$loc.end))
    
    snv.lr <- NULL
    
    if (!is.null(vcf.file)) {
        ov <- findOverlaps(vcf, seg.gr, select = "first")
        
        # Add segment log-ratio to off-target snvs.  For on-target, use observed
        # log-ratio.
        sd.ar <- sd(unlist(geno(vcf)$FA[, tumor.id.in.vcf]))
        snv.lr <- seg$seg.mean[ov]
        ov.vcfexon <- findOverlaps(vcf, exon.gr)
        snv.lr[queryHits(ov.vcfexon)] <- log.ratio[subjectHits(ov.vcfexon)]
        if (sum(is.na(snv.lr)) > 0) {
            n.vcf.before.filter <- nrow(vcf)
            vcf <- vcf[!is.na(snv.lr)]
            mapping.bias <- mapping.bias[!is.na(snv.lr)]
            prior.somatic <- prior.somatic[!is.na(snv.lr)]
            
            # make sure all SNVs are in covered segments
            flog.info("Removing %i variants outside segments.", n.vcf.before.filter - nrow(vcf))
        }
        ov <- findOverlaps(seg.gr, vcf)
        flog.info("Using %i variants.", nrow(vcf))
    }
    
    # get target log-ratios for all segments
    ov.se <- findOverlaps(seg.gr, exon.gr)
    exon.lrs <- lapply(seq_len(nrow(seg)), function(i) log.ratio[subjectHits(ov.se)[queryHits(ov.se) == 
        i]])
    exon.lrs <- lapply(exon.lrs, function(x) subset(x, !is.na(x) & !is.infinite(x)))
    
    # estimate stand. dev. for target logR within targets. this will be used as proxy
    # for sample error.
    targetsPerSegment <- sapply(exon.lrs, length)
    
    if (!sum(targetsPerSegment > 100, na.rm = TRUE)) {
        .stopRuntimeError("Only tiny segments.")
    }
    
    sd.seg <- median(sapply(exon.lrs, sd), na.rm = TRUE)
    
    # if user provided seg file, then we do not have access to the log-ratios and
    # need to use the user provided noise estimate also, don't do outlier smoothing
    # when we use already segmented data
    if (!is.null(seg.file)) {
        sd.seg <- seg.file.sdev
    } else {
    #    exon.lrs <- lapply(exon.lrs, .smoothOutliers)
    }
    
    # renormalize, in case segmentation function changed means
    exon.lrs <- lapply(seq_along(exon.lrs), function(i) {
        if (length(exon.lrs[[i]]) < 3) 
            return(exon.lrs[[i]])
        exon.lrs[[i]] - mean(exon.lrs[[i]]) + seg$seg.mean[i]
    })
    
    max.exon.ratio <- 7
    
    # show log-ratio histogram
    if (plot.cnv) {
        if (!is.null(seg.file)) {
            seg.orig <- read.delim(seg.file)
            par(mfrow = c(2, 1))
            hist(do.call(c, lapply(seq_len(nrow(seg.orig)), function(i) rep(seg.orig$seg.mean[i], 
                seg.orig$num.mark[i]))), breaks = 100, xlab = "log2 ratio", main = paste(sampleid, 
                "(original segmentation)"))
        }
        hist(do.call(c, lapply(seq_len(nrow(seg)), function(i) rep(seg$seg.mean[i], 
            seg$num.mark[i]))), breaks = 100, xlab = "log2 ratio", main = sampleid)
        par(mfrow = c(1, 1))
    }
    
    # initialize variables
    llik <- -Inf
    li <- .getSegSizes(seg)
    C <- rep(2, length(li))
    
    if (sum(li < 0) > 0) 
        .stopRuntimeError("Some segments have negative size.")
    flog.info("Mean standard deviation of log-ratios: %.2f", sd.seg)
    log.ratio.offset <- rep(0, nrow(seg))
    
    flog.info("2D-grid search of purity and ploidy...")
    
    # find local maxima. use a coarser grid for purity, otherwise we will get far too
    # many solutions, which we will need to cluster later anyways.
    if (!is.null(candidates)) {
        candidate.solutions <- candidates
    } else {
        candidate.solutions <- .optimizeGrid(test.purity = seq(max(0.1, min(test.purity)), 
            min(0.99, max(test.purity)), by = 1/30), min.ploidy, max.ploidy, test.num.copy = test.num.copy, 
            exon.lrs, seg, sd.seg, li, max.exon.ratio, max.non.clonal)
        
        # if we have > 20 somatic mutations, we can try estimating purity based on
        # allelic fractions and assuming diploid genomes.
        if (!is.null(vcf.file) && sum(prior.somatic > 0.5, na.rm = TRUE) > 20) {
            somatic.purity <- min(max(test.purity), .calcPuritySomaticVariants(vcf, 
                prior.somatic, tumor.id.in.vcf))
            
            candidate.solutions$candidates <- rbind(candidate.solutions$candidates, 
                c(2, somatic.purity, NA, 2))
        }
    }
    
    if (nrow(candidate.solutions$candidates) > max.candidate.solutions) {
        # test the best solutions and everything close to diploid
        idx.keep <- unique(c(seq_len(max.candidate.solutions), which((candidate.solutions$candidates$tumor.ploidy > 
            1.5 & candidate.solutions$candidates$tumor.ploidy < 2.6))))
        candidate.solutions$candidates <- candidate.solutions$candidates[idx.keep, 
            ]
    }
    flog.info(paste(strwrap(paste("Local optima:\n", paste(round(candidate.solutions$candidates$purity, 
        digits = 2), round(candidate.solutions$candidates$ploidy, digits = 2), 
        sep = "/", collapse = ", "))), collapse = "\n"))
    
    simulated.annealing <- TRUE
    
    .optimizeSolution <- function(cpi) {
        
        max.attempts <- 4
        attempt = 0
        while (attempt < max.attempts) {
            attempt <- attempt + 1
            total.ploidy <- candidate.solutions$candidates$ploidy[cpi]
            p <- candidate.solutions$candidates$purity[cpi]
            
            flog.info("Testing local optimum %i/%i at purity %.2f and total ploidy %.2f...", 
                cpi, nrow(candidate.solutions$candidates), p, total.ploidy)
            
            subclonal <- rep(FALSE, nrow(seg))
            old.llik <- -1
            cnt.llik.equal <- 0
            C.likelihood <- matrix(ncol = length(test.num.copy) + 1, nrow = nrow(seg))
            colnames(C.likelihood) <- c(test.num.copy, "Subclonal")
            for (iter in seq_len(iterations)) {
                # test for convergence
                if (abs(old.llik - llik) < 0.0001) {
                  cnt.llik.equal <- cnt.llik.equal + 1
                }
                
                old.llik <- llik
                if (cnt.llik.equal > 3) 
                  break
                subclonal.f <- length(unlist(exon.lrs[subclonal]))/length(unlist(exon.lrs))
                # should not happen, but sometimes does for very unlikely local optima.
                if (subclonal.f > max.non.clonal + 0.1) 
                  break
                if (iter == 1) 
                  log.ratio.offset <- .sampleOffsetFast(test.num.copy, seg, exon.lrs, 
                    sd.seg, p, C, total.ploidy, max.exon.ratio, simulated.annealing, 
                    log.ratio.calibration)
                
                # in the first iteration, we do not have integer copy numbers yet (corresponding
                # to local optima purity/ploidy)
                if (iter > 1) {
                  # calculate posterior probabilities of all requested purities
                  total.ploidy <- p * (sum(li * (C)))/sum(li) + (1 - p) * 2  #ploidy
                  px.rij <- lapply(test.purity, function(px) vapply(which(!is.na(C)), 
                    function(i) .calcLlikSegment(subclonal = subclonal[i], lr = exon.lrs[[i]] + 
                      log.ratio.offset[i], sd.seg = sd.seg, p = px, Ci = C[i], total.ploidy = total.ploidy, 
                      max.exon.ratio = max.exon.ratio), double(1)))
                  px.rij.s <- sapply(px.rij, sum, na.rm = TRUE) + log(prior.purity)
                  
                  if (simulated.annealing) 
                    px.rij.s <- px.rij.s * exp(iter/4)
                  
                  px.rij.s <- exp(px.rij.s - max(px.rij.s))
                  # Gibbs sample purity
                  p <- test.purity[min(which(runif(n = 1, min = 0, max = sum(px.rij.s)) <= 
                    cumsum(px.rij.s)))]
                  total.ploidy <- p * (sum(li * (C)))/sum(li) + (1 - p) * 2
                  # Gibbs sample offset
                  if (iter > 2) 
                    log.ratio.offset <- .sampleOffset(subclonal, seg, exon.lrs, sd.seg, 
                      p, C, total.ploidy, max.exon.ratio, simulated.annealing, iter, 
                      log.ratio.calibration)
                }
                
                if (!sum(!is.na(C))) {
                  .stopRuntimeError("Could not assign integer copy numbers", " to segments.")
                }
                
                # calculate the log-liklihood of purity and integer copy numbers plus clonal vs
                # subclonal status
                llik <- sum(vapply(which(!is.na(C)), function(i) .calcLlikSegment(subclonal = subclonal[i], 
                  lr = exon.lrs[[i]] + log.ratio.offset[i], sd.seg = sd.seg, p = p, 
                  Ci = C[i], total.ploidy = total.ploidy, max.exon.ratio = max.exon.ratio), 
                  double(1)))
                
                if (is.na(llik)) {
                  .stopRuntimeError("Could not calculate copy number ", "log-likelihood for purity ", 
                    p, " and total ploidy ", total.ploidy, ".")
                }
                
                for (i in seq_len(nrow(seg))) {
                  # Gibbs sample copy number Step 1: calculate log-likelihoods of fits In the first
                  # iteration, we do not have the integer copy numbers yet, so calculate ploidy
                  # only when we have it next time. Now, use the ploidy from the candidate
                  # solution.
                  if (iter > 1) 
                    total.ploidy <- p * (sum(li * (C)))/sum(li) + (1 - p) * 2
                  
                  p.rij <- vapply(test.num.copy, function(Ci) .calcLlikSegment(subclonal = FALSE, 
                    lr = exon.lrs[[i]] + log.ratio.offset[i], sd.seg = sd.seg, p = p, 
                    Ci = Ci, total.ploidy = total.ploidy, max.exon.ratio = max.exon.ratio), 
                    double(1))
                  
                  # calculate tumor ploidy for all possible copy numbers in this segment
                  ploidy <- vapply(test.num.copy, function(Ci) (sum(li[-i] * (C[-i])) + 
                    li[i] * Ci)/sum(li), double(1))
                  
                  # set probability to zero if ploidy is not within requested range
                  log.prior.ploidy <- log(ifelse(ploidy < min.ploidy | ploidy > max.ploidy, 
                    0, 1))
                  if (iter > 1) 
                    p.rij <- p.rij + log.prior.ploidy
                  
                  frac.homozygous.loss <- vapply(test.num.copy, function(Ci) (sum(li[-i] * 
                    ifelse(C[-i] == 0, 1, 0)) + li[i] * ifelse(Ci == 0, 1, 0))/sum(li), 
                    double(1))

                  if (li[i] > max.homozygous.loss[2] && test.num.copy[1]<1) {
                       frac.homozygous.loss[1] <- 1
                  }     
                  log.prior.homozygous.loss <- log(ifelse(frac.homozygous.loss > 
                    max.homozygous.loss[1], 0, 1))
                  if (iter > 1) 
                    p.rij <- p.rij + log.prior.homozygous.loss
                  
                  # model sub clonal state with a uniform distribution
                  p.rij <- c(p.rij, .calcLlikSegmentSubClonal(exon.lrs[[i]] + log.ratio.offset[i], 
                    max.exon.ratio))
                  
                  C.likelihood[i, ] <- exp(p.rij - max(p.rij))
                  
                  if (simulated.annealing) 
                    p.rij <- p.rij * exp(iter/4)
                  
                  # Because we are in log-space, sample relative to most likely fit
                  p.rij <- exp(p.rij - max(p.rij))
                  # Now Gibbs sample best fit
                  z <- runif(n = 1, min = 0, max = sum(p.rij))
                  if (is.na(z)) {
                    message(paste("Iter:", iter, "i:", i, "ploidy:", paste(ploidy, 
                      collapse = "/"), "p.rij:", paste(p.rij)))
                    .stopRuntimeError("Could not fit SCNAs.")
                  }
                  id <- min(which(z <= cumsum(p.rij)))
                  old.C <- C[i]
                  opt.C <- (2^(seg$seg.mean + log.ratio.offset) * total.ploidy)/p - 
                    ((2 * (1 - p))/p)
                  opt.C[opt.C < 0] <- 0
                  if (id > length(test.num.copy)) {
                    # optimal non-integer copy number
                    C[i] <- opt.C[i]
                    subclonal[i] <- TRUE
                  } else {
                    C[i] <- test.num.copy[id]
                    subclonal[i] <- FALSE
                  }
                }
            }
            if (subclonal.f < max.non.clonal && abs(total.ploidy - candidate.solutions$candidates$ploidy[cpi]) < 
                1) 
                break
            log.ratio.calibration <- log.ratio.calibration + 0.25
            if (attempt < max.attempts) {
                flog.info("Recalibrating log-ratios...")
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
                SNV.posterior <- list(beta.model = list(llik = -Inf))
            }
            return(list(log.likelihood = llik, purity = p, ploidy = weighted.mean(C, 
                li), total.ploidy = total.ploidy, seg = seg.adjusted, C.likelihood = data.frame(C.likelihood/rowSums(C.likelihood), 
                ML.C = C, ML.Subclonal = subclonal), SNV.posterior = SNV.posterior, 
                fraction.subclonal = subclonal.f, fraction.homozygous.loss = sum(li[which(C < 
                  0.01)])/sum(li), gene.calls = NA, log.ratio.offset = log.ratio.offset, 
                SA.iterations = iter))
        }
        
        if (!is.null(gc.gene.file)) {
            gene.calls <- .getGeneCalls(seg.adjusted, gc.data, log.ratio, fun.focal, 
                args.focal, chr.hash)
        } else {
            gene.calls <- NA
        }
        
        if (!is.null(vcf.file)) {
            if (post.optimize) {
                idx <- c(match(p, test.purity), 
                    (max(1, match(p, test.purity) - 4)):(min(length(test.purity), 
                     match(p, test.purity) + 4)))
                idx <- idx[!duplicated(idx)]
                tp <- test.purity[idx]
                pp <- prior.purity[idx]
            } else {
                tp <- p
                pp <- 1
            }
            cont.rate <- prior.contamination
            .fitSNV <- function(tp, pp) {
                .fitSNVp <- function(px, cont.rate=prior.contamination) {
                    flog.info("Fitting SNVs for purity %.2f, tumor ploidy %.2f and contamination %.2f.", 
                        px, weighted.mean(C, li), cont.rate)
                  
                  list(beta.model = .calcSNVLLik(vcf, tumor.id.in.vcf, ov, px, test.num.copy, 
                    C.likelihood, C, opt.C, snv.model = "beta", prior.somatic, mapping.bias, 
                    snv.lr, sampleid, cont.rate = cont.rate, prior.K = prior.K, 
                    max.coverage.vcf = max.coverage.vcf, non.clonal.M = non.clonal.M, 
                    model.homozygous = model.homozygous, error = error, max.mapping.bias = max.mapping.bias))
                }

                res.snvllik <- lapply(tp[1], .fitSNVp)
                
                if (post.optimize && length(tp)>1) {
                    GoF <- .getGoF(list(SNV.posterior=res.snvllik[[1]]))
                    idx <- 1
                    if (GoF < (min.gof-0.05)) { 
                        flog.info("Poor goodness-of-fit (%.3f). Skipping post-optimization.", GoF)
                    } else {    
                        res.snvllik <- c(res.snvllik, lapply(tp[-1], .fitSNVp))
                      px.rij <- lapply(tp, function(px) vapply(which(!is.na(C)), function(i) .calcLlikSegment(subclonal = subclonal[i], 
                        lr = exon.lrs[[i]] + log.ratio.offset[i], sd.seg = sd.seg, p = px, 
                        Ci = C[i], total.ploidy = px * (sum(li * (C)))/sum(li) + (1 - 
                          px) * 2, max.exon.ratio = max.exon.ratio), double(1)))
                      
                      px.rij.s <- sapply(px.rij, sum, na.rm = TRUE) + log(pp) + vapply(res.snvllik, 
                        function(x) x$beta.model$llik, double(1))
                      idx <- which.max(px.rij.s)
                  }
                } else {
                  idx <- 1
                }
                p <- tp[idx]
                flog.info("Optimized purity: %.2f", p)
                SNV.posterior <- res.snvllik[[idx]]
                list(p = p, SNV.posterior = SNV.posterior, llik = px.rij.s[idx])
            }
            fitRes <- .fitSNV(tp, pp)
            p <- fitRes$p
            SNV.posterior <- fitRes$SNV.posterior
        }
         
        list(log.likelihood = llik, purity = p, ploidy = weighted.mean(C, li),
            total.ploidy = total.ploidy, seg = seg.adjusted, 
            C.posterior = data.frame(C.likelihood/rowSums(C.likelihood), ML.C =
            C, Opt.C = opt.C, ML.Subclonal = subclonal),
            C.likelihood=C.likelihood, SNV.posterior = SNV.posterior,
            fraction.subclonal = subclonal.f, fraction.homozygous.loss =
            sum(li[which(C < 0.01)])/sum(li), gene.calls = gene.calls,
            log.ratio.offset = log.ratio.offset, SA.iterations = iter, failed =
            FALSE)
    }
    
    results <- lapply(seq_len(nrow(candidate.solutions$candidates)), .optimizeSolution)
    
    results <- .rankResults(results)
    results <- .filterDuplicatedResults(results)
    cont.rate <- prior.contamination

    if (grepl("CONTAMINATION", vcf.filtering$flag_comment)) {
        cont.rate <- .plotContamination(
            results[[1]]$SNV.posterior$beta.model$posteriors, 
            max.mapping.bias, plot=FALSE)
        if (cont.rate > prior.contamination) {
            flog.info("Initial guess of contamination rate: %.3f", cont.rate)
        }    
    }
    ## optimize contamination. we just re-run the fitting 
    if (grepl("CONTAMINATION", vcf.filtering$flag_comment) && 
        cont.rate>prior.contamination) {
        flog.info("Optimizing contamination rate...")
             
        res.snvllik <-
            .calcSNVLLik(vcf, tumor.id.in.vcf,
                  ov, results[[1]]$purity, test.num.copy, results[[1]]$C.likelihood, 
                  results[[1]]$C.posterior$ML.C,
                  results[[1]]$C.posterior$Opt.C,
                  snv.model = "beta", prior.somatic, mapping.bias,
                  snv.lr, sampleid, cont.rate = cont.rate, prior.K = prior.K,
                  max.coverage.vcf = max.coverage.vcf, non.clonal.M = non.clonal.M,
                  model.homozygous = model.homozygous, error = error,
                  max.mapping.bias = max.mapping.bias)
        results[[1]]$SNV.posterior$beta.model <- res.snvllik
        cont.rate <- .plotContamination(
                    results[[1]]$SNV.posterior$beta.model$posteriors,
                                max.mapping.bias, plot=FALSE)
                flog.info("Optimized contamination rate: %.3f", cont.rate)
        results[[1]]$SNV.posterior$beta.model$posterior.contamination <- cont.rate
        # undo flagging if contamination rate is low
    }
    if (cont.rate<0.001 && vcf.filtering$flag_comment == 
        "POTENTIAL SAMPLE CONTAMINATION") {
        vcf.filtering$flag_comment <- ""
        vcf.filtering$flag <- FALSE
    }    
    results <- .flagResults(results, max.non.clonal = max.non.clonal, max.logr.sdev = max.logr.sdev, 
        logr.sdev = sd.seg, max.segments = max.segments, min.gof = min.gof, flag = vcf.filtering$flag, 
        flag_comment = vcf.filtering$flag_comment, dropout = dropoutWarning, use.somatic.status = args.filterVcf$use.somatic.status, 
        model.homozygous = model.homozygous)
    
    if (!is.null(normalDB) && length(results) && !is.null(results[[1]]$gene.calls)) {
        results <- .addVoomToGeneCalls(results, tumor.coverage.file, normalDB, 
            gc.gene.file)
    }
    
    if (length(results) < 1) {
        flog.warn("Could not find valid purity and ploidy solution.")
    }
    .logFooter()
    list(candidates = candidate.solutions, results = results, input = list(tumor = tumor.coverage.file, 
        normal = normal.coverage.file, log.ratio = data.frame(probe = normal[, 1], 
            log.ratio = log.ratio), log.ratio.sdev = sd.seg, vcf = vcf, sampleid = sampleid, 
        sex = sex, sex.vcf = sex.vcf, chr.hash = chr.hash, centromeres = centromeres))
}


