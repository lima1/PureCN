#' Run PureCN implementation of ABSOLUTE
#'
#' This function takes as input tumor and normal control coverage data and
#' a VCF containing allelic fractions of germline variants and somatic
#' mutations. Normal control does not need to be from the same patient.
#' In case VCF does not contain somatic status, it should contain dbSNP and
#' optionally COSMIC annotation. Returns purity and ploidy combinations,
#' sorted by likelihood score. Provides copy number and LOH data, by both
#' gene and genomic region.
#'
#'
#' @param normal.coverage.file Coverage file of normal control (optional
#' if log.ratio is provided - then it will be only used to filter low coverage
#' exons).  Should be already GC-normalized with
#' \code{\link{correctCoverageBias}}.  Needs to be either a file name or data
#' read with the \code{\link{readCoverageFile}} function.
#' @param tumor.coverage.file Coverage file of tumor. If \code{NULL},
#' requires \code{seg.file} and an interval file via \code{interval.file}.
#' Should be already GC-normalized with \code{\link{correctCoverageBias}}.
#' Needs to be either a file name or data read with the
#' \code{\link{readCoverageFile}} function.
#' @param log.ratio Copy number log-ratios for all exons in the coverage files.
#' If \code{NULL}, calculated based on coverage files.
#' @param seg.file Segmented data. Optional, to support third-pary
#' segmentation tools.  If \code{NULL}, use coverage files or
#' \code{log.ratio} to segment the data.
#' @param seg.file.sdev If \code{seg.file} provided, the log-ratio standard
#' deviation, used to model likelihood of sub-clonal copy number events.
#' @param vcf.file VCF file.
#' Optional, but typically needed to select between local optima of similar
#' likelihood. Can also be a \code{CollapsedVCF}, read with the \code{readVcf}
#' function.  Requires a DB info flag for dbSNP membership. The default
#' \code{fun.setPriorVcf} function will also look for a Cosmic.CNT slot (see
#' \code{cosmic.vcf.file}), containing the hits in the COSMIC database. Again,
#' do not expect very useful results without a VCF file.
#' @param normalDB Normal database, created with
#' \code{\link{createNormalDatabase}}. If provided, used to calculate gene-level
#' p-values (requires \code{Gene} column in \code{interval.file}) and to filter
#' targets with low coverage in the pool of normal samples.
#' @param genome Genome version, for example hg19. See \code{readVcf}.
#' @param centromeres A \code{GRanges} object with centromere positions.
#' If \code{NULL}, use pre-stored positions for genome versions
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
#' @param fun.filterIntervals Function for filtering low-quality intervals in the
#' coverage files. Needs to return a \code{logical} vector whether an interval
#' should be used for segmentation. Defaults to \code{\link{filterIntervals}}.
#' @param args.filterIntervals Arguments for target filtering function. Arguments
#' \code{normal}, \code{tumor}, \code{log.ratio},
#' \code{min.coverage}\code{seg.file} and \code{normalDB} are required and
#' automatically set.
#' @param fun.segmentation Function for segmenting the copy number log-ratios.
#' Expected return value is a \code{data.frame} representation of the
#' segmentation. Defaults to \code{\link{segmentationCBS}}.
#' @param args.segmentation Arguments for segmentation function. Arguments
#' \code{normal}, \code{tumor}, \code{log.ratio}, \code{plot.cnv},
#' \code{sampleid}, \code{vcf}, \code{tumor.id.in.vcf},
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
#' @param min.coverage Minimum coverage in both normal and tumor. Intervals and
#' variants with lower coverage are ignored. This value is provided to the
#' \code{args.filterIntervals} and \code{args.filterVcf} lists, but can be
#' overwritten in these lists if different cutoffs for the coverage and variant
#' filters are wanted. To increase the sensitivity of homozygous deletions in
#' high purity samples, the coverage cutoff in tumor is automatically lowered by
#' 50 percent if the normal coverage is high.
#' @param max.coverage.vcf This will set the maximum number of reads in the SNV
#' fitting.  This is to avoid that small non-reference biases that come
#' apparent only at high coverages have a dramatic influence on likelihood
#' scores. Only relevant for \code{model = "beta"}.
#' @param max.non.clonal Maximum genomic fraction assigned to a subclonal copy
#' number state.
#' @param max.homozygous.loss \code{double(2)} with maximum chromosome fraction
#' assigned to homozygous loss and maximum size of a homozygous loss segment.
#' @param non.clonal.M Average expected cellular fraction of sub-clonal somatic
#' mutations. This is to calculate expected allelic fractions of a single
#' sub-clonal bin for variants. For all somatic variants, more accurate cellular
#' fractions are calculated.
#' @param max.mapping.bias Exclude variants with high mapping bias from the
#' likelihood score calculation. Note that bias is reported on an inverse
#' scale; a variant with mapping bias of 1 has no bias.
#' @param max.pon Exclude variants found more than \code{max.pon} times in
#' pool of normals and not in dbSNP. Requires \code{mapping.bias.file} in
#' \code{\link{setMappingBiasVcf}}. Should be set to a value high enough
#' to be much more likely an artifact and not a true germline variant not
#' present in dbSNP.
#' @param min.variants.segment Flag segments with fewer variants. The
#' minor copy number estimation is not reliable with insufficient variants.
#' @param iterations Maximum number of iterations in the Simulated Annealing
#' copy number fit optimization. Note that this an integer optimization problem
#' that should converge quickly. Allowed range is 10 to 250.
#' @param log.ratio.calibration Re-calibrate log-ratios in the window
#' \code{purity*log.ratio.calibration}.
#' @param smooth.log.ratio Smooth \code{log.ratio} using the \code{DNAcopy}
#' package.
#' @param model.homozygous Homozygous germline SNPs are uninformative and by
#' default removed. In 100 percent pure samples such as cell lines, however,
#' heterozygous germline SNPs appear homozygous in case of LOH. Setting this
#' parameter to \code{TRUE} will keep homozygous SNPs and include a homozygous
#' SNP state in the likelihood model. Not necessary when matched normal samples
#' are available.
#' @param error Estimated sequencing error rate. Used to calculate minimum
#' number of supporting reads for variants using
#' \code{\link{calculatePowerDetectSomatic}}. Also used to calculate the
#' probability of homozygous SNP allelic fractions (assuming reference reads
#' are sequencing errors).
#' @param interval.file A mapping file that assigns GC content and gene symbols
#' to each exon in the coverage files. Used for generating gene-level calls.
#' First column in format CHR:START-END. Second column GC content (0 to 1).
#' Third column gene symbol. This file is generated with the
#' \code{\link{preprocessIntervals}} function.
#' @param max.dropout Measures GC bias as ratio of coverage in AT-rich (GC <
#' 0.5) versus GC-rich on-target regions (GC >= 0.5). High coverage drop-out might
#' indicate that  data was not GC-normalized (optional with larger pool of
#' normal samples). A warning pointing to a normalized log-ratio drop-out likely
#' indicates that the sample quality is insufficient. For log-ratio drop-out,
#' a warning is thrown when half the \code{max.dropout} is reached since it
#' is calculated using both tumor and normal.
#' Requires \code{interval.file}.
#' @param min.logr.sdev Minimum log-ratio standard deviation used in the
#' model. Useful to make fitting more robust to outliers in very clean
#' data.
#' @param max.logr.sdev Flag noisy samples with segment log-ratio standard
#' deviation larger than this. Assay specific and needs to be calibrated.
#' @param max.segments Flag noisy samples with a large number of segments.
#' Assay specific and needs to be calibrated.
#' @param min.gof Flag purity/ploidy solutions with poor fit.
#' @param min.variants Do not attempt to fit allelic fractions for samples
#' with fewer variants passing all filters.
#' @param plot.cnv Generate segmentation plots.
#' @param vcf.field.prefix Prefix all newly created VCF field names with
#' this string.
#' @param cosmic.vcf.file Add a \code{Cosmic.CNT} info field to the provided
#' \code{vcf.file} using a VCF file containing the COSMIC database. The default
#' \code{fun.setPriorVcf} function will give variants found in the COSMIC database
#' a higher prior probability of being somatic. Not used in likelhood model
#' when matched normal is available in \code{vcf.file}. Should be compressed
#' and indexed with bgzip and tabix, respectively.
#' @param DB.info.flag Flag in INFO of VCF that marks presence in common
#' germline databases. Defaults to \code{DB} that may contain somatic variants
#' if it is from an unfiltered dbSNP VCF.
#' @param POPAF.info.field As alternative to a flag, use an info field that
#' contains population allele frequencies. The \code{DB} info flag has priority
#' over this field when both exist.
#' @param Cosmic.CNT.info.field Info field containing hits in the Cosmic database
#' @param min.pop.af Minimum population allele frequency in
#' \code{POPAF.info.field} to set a high germline prior probability.
#' @param model Use either a beta or a beta-binomial distribution for fitting
#' observed to expected allelic fractions of alterations in \code{vcf.file}.
#' The latter can be useful to account for significant overdispersion, for example
#' due to mapping biases when no pool of normals is available or due to other
#' unmodeled biases, e.g. amplification biases. The beta-binomial model is only
#' recommended with a sufficiently sized pool of normal samples
#' (more than 10 normals)
#' @param post.optimize Optimize purity using final SCNA-fit and variants. This
#' might take a long time when lots of variants need to be fitted, but will
#' typically result in a slightly more accurate purity, especially for rather
#' silent genomes or very low purities. Otherwise, it will just use the purity
#' determined via the SCNA-fit.
#' @param speedup.heuristics Tries to avoid spending computation time on
#' local optima that are unlikely correct. Set to 0 to turn this off, to 1 to
#' only apply heuristics that in worst case will decrease accuracy slightly or
#' to 2 to turn on all heuristics.
#' @param BPPARAM \code{BiocParallelParam} object. If \code{NULL}, does not
#' use parallelization for fitting local optima.
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
#' normal.coverage.file <- system.file('extdata', 'example_normal_tiny.txt',
#'     package = 'PureCN')
#' tumor.coverage.file <- system.file('extdata', 'example_tumor_tiny.txt',
#'     package = 'PureCN')
#' vcf.file <- system.file('extdata', 'example.vcf.gz',
#'     package = 'PureCN')
#' interval.file <- system.file('extdata', 'example_intervals_tiny.txt',
#'     package = 'PureCN')
#'
#' # The max.candidate.solutions, max.ploidy and test.purity parameters are set to
#' # non-default values to speed-up this example.  This is not a good idea for real
#' # samples.
#' ret <-runAbsoluteCN(normal.coverage.file = normal.coverage.file,
#'     tumor.coverage.file = tumor.coverage.file, genome = 'hg19',
#'     vcf.file = vcf.file, sampleid = 'Sample1',
#'     interval.file = interval.file, max.ploidy = 4,
#'     test.purity = seq(0.3, 0.7, by = 0.05), max.candidate.solutions = 1)
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
#' res <- runAbsoluteCN(seg.file = seg.file, fun.segmentation = funSeg,
#'     max.ploidy = 4, test.purity = seq(0.3, 0.7, by = 0.05),
#'     max.candidate.solutions = 1,
#'     genome='hg19', interval.file = interval.file)
#'
#' @export runAbsoluteCN
#' @import DNAcopy
#' @import IRanges
#' @import VariantAnnotation
#' @importFrom Biobase package.version
#' @importFrom GenomicRanges GRanges tile
#' @importFrom stats complete.cases dbeta dnorm dunif runif weighted.mean
#'             dbinom C
#' @importFrom utils data read.delim tail packageVersion object.size
#' @importFrom S4Vectors queryHits subjectHits DataFrame expand %in%
#' @importFrom data.table data.table
#' @importFrom futile.logger flog.info flog.warn flog.fatal flog.debug
#'             flog.threshold flog.appender appender.tee
#' @importFrom VGAM dbetabinom.ab
runAbsoluteCN <- function(normal.coverage.file = NULL,
    tumor.coverage.file = NULL, log.ratio = NULL, seg.file = NULL,
    seg.file.sdev = 0.4, vcf.file = NULL, normalDB = NULL, genome,
    centromeres = NULL, sex = c("?", "F", "M", "diploid"),
    fun.filterVcf = filterVcfMuTect, args.filterVcf = list(),
    fun.setPriorVcf = setPriorVcf, args.setPriorVcf = list(),
    fun.setMappingBiasVcf = setMappingBiasVcf, args.setMappingBiasVcf = list(),
    fun.filterIntervals = filterIntervals, args.filterIntervals = list(),
    fun.segmentation = segmentationCBS, args.segmentation = list(),
    fun.focal = findFocal, args.focal = list(),
    sampleid = NULL, min.ploidy = 1.4, max.ploidy = 6, test.num.copy = 0:7,
    test.purity = seq(0.15, 0.95, by = 0.01), prior.purity = NULL,
    prior.K = 0.999, prior.contamination = 0.01, max.candidate.solutions = 20,
    candidates = NULL, min.coverage = 15, max.coverage.vcf = 300,
    max.non.clonal = 0.2, max.homozygous.loss = c(0.05, 1e07),
    non.clonal.M = 1 / 3, max.mapping.bias = 0.8, max.pon = 3, iterations = 30,
    min.variants.segment = 5, log.ratio.calibration = 0.1, smooth.log.ratio = TRUE,
    model.homozygous = FALSE, error = 0.001,
    interval.file = NULL, max.dropout = c(0.95, 1.1),
    min.logr.sdev = 0.15, max.logr.sdev = 0.6,
    max.segments = 300, min.gof = 0.8, min.variants = 20,
    plot.cnv = TRUE,
    vcf.field.prefix = "",
    cosmic.vcf.file = NULL, 
    DB.info.flag = "DB",
    POPAF.info.field = "POP_AF",
    Cosmic.CNT.info.field = "Cosmic.CNT",
    min.pop.af = 0.001,
    model = c("beta", "betabin"),
    post.optimize = FALSE, speedup.heuristics = 2, BPPARAM = NULL,
    log.file = NULL, verbose = TRUE) {

    if (!verbose) flog.threshold("WARN")
    if (!is.null(log.file)) flog.appender(appender.tee(log.file))
    
     # log function arguments
    try(.logHeader(as.list(match.call())[-1]), silent = TRUE)

    model <- match.arg(model)

    # defaults to equal priors for all tested purity values
    if (is.null(prior.purity)) {
        prior.purity <- rep(1, length(test.purity)) / length(test.purity)
    }
    # argument checking
    .checkParameters(test.purity, min.ploidy, max.ploidy, max.non.clonal,
        max.homozygous.loss, sampleid, prior.K, prior.contamination, prior.purity,
        iterations, min.gof, model.homozygous, interval.file, log.ratio.calibration,
        test.num.copy, max.mapping.bias)

    test.num.copy <- sort(test.num.copy)

    if (!is.null(BPPARAM) && !requireNamespace("BiocParallel", quietly = TRUE)) {
        flog.warn("Install BiocParallel for parallel optimization.")
        BPPARAM <- NULL
    } else if (!is.null(BPPARAM)) {
        flog.info("Using BiocParallel for parallel optimization.")
    }

    flog.info("Loading coverage files...")

    if (!is.null(normal.coverage.file)) {
        if (is.character(normal.coverage.file)) {
            normal.coverage.file <- normalizePath(normal.coverage.file, mustWork = TRUE)
            normal <- readCoverageFile(normal.coverage.file)
        } else {
            normal <- normal.coverage.file
        }
    }
    if (is.null(tumor.coverage.file)) {
        if ((is.null(seg.file) && is.null(log.ratio)) ||
            is.null(interval.file)) {
            .stopUserError("Missing tumor.coverage.file requires seg.file or ",
                           "log.ratio and interval.file.")
        }
        tumor <- .gcGeneToCoverage(interval.file, min.coverage + 1)
        tumor.coverage.file <- tumor
    } else if (is.character(tumor.coverage.file)) {
        tumor.coverage.file <- normalizePath(tumor.coverage.file, mustWork = TRUE)
        tumor <- readCoverageFile(tumor.coverage.file)
        if (is.null(sampleid))
            sampleid <- basename(tumor.coverage.file)
    } else {
        tumor <- tumor.coverage.file
        if (is.null(sampleid))
            sampleid <- "Sample.1"
    }

    chr.hash <- .getChrHash(seqlevels(tumor))

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
            if (is.null(normal.coverage.file)) {
                log.ratio <- .createFakeLogRatios(tumor, seg.file, sampleid,
                                                  chr.hash, model.homozygous, max.logr.sdev)
                smooth.log.ratio <- FALSE
            } else {
                flog.info("seg.file and normal.coverage.file provided. Using both.")
                log.ratio <- calculateLogRatio(normal, tumor)
            }
            if (is.null(sampleid))
                sampleid <- read.delim(seg.file, as.is = TRUE)[1, 1]
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
        if (length(log.ratio) != length(tumor)) {
            .stopUserError("Length of log.ratio different from tumor ",
                           "coverage.")
        }
        if (is.null(sampleid))
            sampleid <- "Sample.1"
    }

    sex <- .getSex(match.arg(sex), normal, tumor)
    tumor <- .fixAllosomeCoverage(sex, tumor)

    if (!is.null(interval.file)) {
        tumor <- .addGCData(tumor, interval.file)
    }
    if (is.null(centromeres) && !missing(genome)) {
        centromeres <- .getCentromerePositions(centromeres, genome,
            if (is.null(tumor)) NULL else .getSeqlevelsStyle(tumor))
    }

    args.filterIntervals <- c(list(normal = normal, tumor = tumor,
        log.ratio = log.ratio, seg.file = seg.file,
        normalDB = normalDB), args.filterIntervals)

    # make it possible to provide different coverages for the different
    # filters
    if (is.null(args.filterIntervals$min.coverage)) {
        args.filterIntervals$min.coverage <- min.coverage
    }

    intervalsUsed <- do.call(fun.filterIntervals,
        .checkArgs(args.filterIntervals, "filterIntervals"))

    # chr.hash is an internal data structure, so we need to do this separately.
    intervalsUsed <- .filterIntervalsChrHash(intervalsUsed, tumor, chr.hash)
    intervalsUsed <- .filterIntervalsCentromeres(intervalsUsed, tumor, centromeres)
    intervalsUsed <- which(intervalsUsed)
    if (length(tumor) != length(normal) ||
        length(tumor) != length(log.ratio)) {
        .stopRuntimeError("Mis-aligned coverage data.")
    }
    tumor <- tumor[intervalsUsed, ]
    normal <- normal[intervalsUsed, ]
    log.ratio <- log.ratio[intervalsUsed]
    flog.info("Using %i intervals (%i on-target, %i off-target).", length(tumor),
        sum(tumor$on.target, na.rm = TRUE),
        sum(!tumor$on.target, na.rm = TRUE))

    if (!length(tumor)) {
        .stopUserError("No intervals passing filters.")
    }
    if (!sum(!tumor$on.target, na.rm = TRUE)) {
        flog.info("No off-target intervals. If this is hybrid-capture data,%s",
            " consider adding them.")
    } else {
         flog.info("Ratio of mean on-target vs. off-target read counts: %.2f",
            mean(tumor$counts[tumor$on.target], na.rm = TRUE) /
            mean(tumor$counts[!tumor$on.target], na.rm = TRUE))
         
         flog.info("Mean off-target bin size: %.0f",
            mean(width(tumor[!tumor$on.target]), na.rm = TRUE))
    }
    if (!is.null(normalDB$sd$weights)) {
        tumor$weights <- subsetByOverlaps(normalDB$sd$weights, tumor)$weights
    }
    # not needed anymore
    normalDB <- NULL

    if (smooth.log.ratio) {
        log.ratio <- smooth.CNA(
            .getCNAobject(log.ratio, normal, chr.hash, "sample"))$sample
    }

    dropoutWarning <- FALSE
    # clean up noisy targets, but not if the segmentation was already provided.
    if (is.null(seg.file)) {
        if (!is.null(interval.file)) {
            dropoutWarning <- .checkGCBias(normal, tumor, log.ratio, max.dropout)
        } else {
            flog.info("No interval.file provided. Cannot check for any GC-biases.")
        }
    }

    vcf <- NULL
    vcf.germline <- NULL
    tumor.id.in.vcf <- NULL
    normal.id.in.vcf <- NULL
    vcf.filtering <- list(flag = FALSE, flag_comment = "")
    sex.vcf <- NULL

    if (!is.null(vcf.file)) {
        flog.info("Loading VCF...")
        vcf <- .readAndCheckVcf(vcf.file, genome = genome,
            DB.info.flag = DB.info.flag, POPAF.info.field = POPAF.info.field,
            min.pop.af = min.pop.af, error = error,
            vcf.field.prefix = vcf.field.prefix)

        if (length(intersect(seqlevels(tumor), seqlevels(vcf))) < 1) {
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
            sex.vcf <- getSexFromVcf(vcf, tumor.id.in.vcf,
                use.somatic.status = args.filterVcf$use.somatic.status)
            if (!is.na(sex.vcf) && sex %in% c("F", "M") && sex.vcf != sex) {
                flog.warn("Sex mismatch of coverage and VCF. %s%s",
                        "Could be because of noisy data, contamination, ",
                  "loss of chrY or a mis-alignment of coverage and VCF.")
            }
        }
        n.vcf.before.filter <- nrow(vcf)

        args.filterVcf <- c(list(vcf = vcf, tumor.id.in.vcf = tumor.id.in.vcf,
            model.homozygous = model.homozygous, error = error,
            DB.info.flag = DB.info.flag,
            target.granges = tumor[tumor$on.target]), args.filterVcf)
        if (is.null(args.filterVcf$min.coverage)) {
            args.filterVcf$min.coverage <- min.coverage
        }
        vcf.filtering <- do.call(fun.filterVcf,
            .checkArgs(args.filterVcf, "filterVcf"))

        vcf <- vcf.filtering$vcf

        if (!is.null(cosmic.vcf.file)) {
            vcf <- .addCosmicCNT(vcf, cosmic.vcf.file, Cosmic.CNT.info.field)
        }

        args.setPriorVcf <- c(list(vcf = vcf, tumor.id.in.vcf = tumor.id.in.vcf,
            DB.info.flag = DB.info.flag,
            Cosmic.CNT.info.field = Cosmic.CNT.info.field), args.setPriorVcf)
        vcf <- do.call(fun.setPriorVcf,
            .checkArgs(args.setPriorVcf, "setPriorVcf"))

        # get mapping bias
        args.setMappingBiasVcf$vcf <- vcf
        args.setMappingBiasVcf$tumor.id.in.vcf <- tumor.id.in.vcf
        vcf <- do.call(fun.setMappingBiasVcf,
            .checkArgs(args.setMappingBiasVcf, "setMappingBiasVcf"))
        idxHqGermline <- info(vcf)[[paste0(vcf.field.prefix, "PR")]] < 0.1 &
            info(vcf)[[paste0(vcf.field.prefix, "MBB")]] >= max.mapping.bias & 
            info(vcf)[[paste0(vcf.field.prefix, "MBB")]] <= (2 - max.mapping.bias)

        flog.info("Excluding %i novel or poor quality variants from segmentation.", sum(!idxHqGermline))
        # for larger pool of normals, require that we have seen the SNP
        pon.count <- info(vcf)[[paste0(vcf.field.prefix, "MBPON")]]
        if (any(!is.na(pon.count)) &&
            max(pon.count, na.rm = TRUE) > 10) {
            idxHqGermline <- idxHqGermline & pon.count > 0
            flog.info("Excluding %i variants not in pool of normals from segmentation.",
                 sum(!pon.count > 0))
        }
        vcf.germline <- vcf[idxHqGermline]
    }

    flog.info("Sample sex: %s", sex)
    flog.info("Segmenting data...")

    segProvided <- readSegmentationFile(seg.file, sampleid, model.homozygous = model.homozygous)
    
    args.segmentation <- c(list(normal = normal, tumor = tumor, log.ratio = log.ratio, 
        seg = segProvided, plot.cnv = plot.cnv,
        sampleid = sampleid, vcf = vcf.germline, tumor.id.in.vcf = tumor.id.in.vcf,
        normal.id.in.vcf = normal.id.in.vcf, max.segments = max.segments, chr.hash = chr.hash,
        min.logr.sdev = min.logr.sdev, centromeres = centromeres), args.segmentation)

    vcf.germline <- NULL
    seg <- do.call(fun.segmentation,
            .checkArgs(args.segmentation, "segmentation"))

    if (!is.null(seg.file)) {
        seg <- .fixAllosomeSegmentation(sex, seg)
    }

    if (speedup.heuristics > 1) {
        ds <- .getSizeDomState(seg)
        if (ds$fraction.genome > 0.5 && abs(ds$seg.mean) < 0.1) {
            min.ploidy <- 1.5
            max.ploidy <- 3
            flog.info("Highly dominant copy number state with small log-ratio.%s",
                " Skipping search for high ploidy solutions.")
        }
    }
    seg.gr <- GRanges(seqnames = .add.chr.name(seg$chrom, chr.hash),
                IRanges(start = round(seg$loc.start), end = seg$loc.end))

    missing.sl <- setdiff(seqlevelsInUse(tumor), seqlevelsInUse(seg.gr))
    missing.sl <- missing.sl[!missing.sl %in% .getSexChr(seqlevelsInUse(tumor))]

    if (!is.null(seg.file) && length(missing.sl)) {
        .stopUserError("Seqlevels missing in provided segmentation: ",
                       paste(missing.sl, collapse = ","))
    }    

    flog.info("Found %i segments with median size of %.2fMb.",
              length(seg.gr), median(width(seg.gr) / 1e+6))

    snv.lr <- NULL

    if (!is.null(vcf.file)) {
        ov <- findOverlaps(vcf, seg.gr, select = "first")

        # Add segment log-ratio to off-target snvs.  For on-target, use observed
        # log-ratio.
        snv.lr <- seg$seg.mean[ov]
        ov.vcfexon <- findOverlaps(vcf, tumor)
        snv.lr[queryHits(ov.vcfexon)] <- log.ratio[subjectHits(ov.vcfexon)]
        if (anyNA(snv.lr)) {
            n.vcf.before.filter <- .countVariants(vcf)
            vcf <- .removeVariants(vcf, is.na(snv.lr), "segmentation")
            # make sure all variants are in covered segments
            flog.info("Removing %i variants outside segments.", n.vcf.before.filter - .countVariants(vcf))
        }
        ov <- findOverlaps(seg.gr, vcf)
        flog.info("Using %i variants.", .countVariants(vcf))
    }

    # get target log-ratios for all segments
    exon.lrs <- .getExonLrs(seg.gr, tumor, log.ratio)

    sd.seg <- max(median(vapply(exon.lrs, sd, double(1)), na.rm = TRUE),
                  min.logr.sdev)

    if (sum(!tumor$on.target, na.rm = TRUE)) {
        sd.seg.ontarget <- median(vapply(
            .getExonLrs(seg.gr, tumor, log.ratio, idx = tumor$on.target),
            sd, double(1)), na.rm = TRUE)

        sd.seg.offtarget <- median(vapply(
            .getExonLrs(seg.gr, tumor, log.ratio, idx = !tumor$on.target),
            sd, double(1)), na.rm = TRUE)
    }
    # if user provided seg file, then we do not have access to the log-ratios and
    # need to use the user provided noise estimate also, don't do outlier smoothing
    # when we use already segmented data
    if (!is.null(seg.file) && .isFakeLogRatio(log.ratio)) {
        sd.seg <- seg.file.sdev
    }

    # renormalize, in case segmentation function changed means
    exon.lrs <- .postprocessLogRatios(exon.lrs, seg$seg.mean)

    max.exon.ratio <- 4

    # show log-ratio histogram and on/off-target log-ratios
    if (plot.cnv) {
        if (!is.null(segProvided)) {
            par(mfrow = c(2, 1))
            hist(do.call(c, lapply(seq_len(nrow(segProvided)), function(i) rep(segProvided$seg.mean[i],
                segProvided$num.mark[i]))), breaks = 100, xlab = "log2 ratio", main = paste(sampleid,
                "(original segmentation)"))
        }
        hist(do.call(c, lapply(seq_len(nrow(seg)), function(i) rep(seg$seg.mean[i],
            seg$num.mark[i]))), breaks = 100, xlab = "log2 ratio", main = sampleid)
        par(mfrow = c(1, 1))
        .plotLogRatios(log.ratio, tumor$on.target)
    }

    # initialize variables
    llik <- -Inf
    li <- .getSegSizes(seg)
    C <- rep(2, length(li))
    mapd <- list(all = median(abs(diff(log.ratio))))
    if (sum(li < 0) > 0)
        .stopRuntimeError("Some segments have negative size.")
    flog.info("Mean standard deviation of log-ratios: %.2f (MAPD: %.2f)",
        sd.seg, mapd$all)
    if (sum(!tumor$on.target, na.rm = TRUE)) {
        mapd$on.target <- median(abs(diff(log.ratio[tumor$on.target])))
        mapd$off.target <- median(abs(diff(log.ratio[!tumor$on.target])))
        flog.info("Mean standard deviation of on-target log-ratios only: %.2f (MAPD: %.2f)",
            sd.seg.ontarget, mapd$on.target)
        flog.info("Mean standard deviation of off-target log-ratios only: %.2f (MAPD: %.2f)",
            sd.seg.offtarget, mapd$off.target)
    }
    log.ratio.offset <- rep(0, nrow(seg))

    flog.info("2D-grid search of purity and ploidy...")

    # find local maxima. use a coarser grid for purity, otherwise we will get far too
    # many solutions, which we will need to cluster later anyways.
    if (!is.null(candidates)) {
        candidate.solutions <- candidates
    } else {
        candidate.solutions <- .optimizeGrid(.get2DPurityGrid(test.purity),
            min.ploidy, max.ploidy, test.num.copy = test.num.copy,
            exon.lrs, seg, sd.seg, li, max.exon.ratio, max.non.clonal, BPPARAM)

        # if we have > 20 somatic mutations, we can try estimating purity based on
        # allelic fractions and assuming diploid genomes.
        if (!is.null(vcf.file) &&
            sum(info(vcf)[[paste0(vcf.field.prefix, "PR")]] > 0.5,
                na.rm = TRUE) > 20) {
            somatic.purity <- .calcPuritySomaticVariants(vcf, tumor.id.in.vcf)
            somatic.purity <- min(max(test.purity), somatic.purity)
            somatic.purity <- max(min(test.purity), somatic.purity)

            candidate.solutions$candidates <- .filterDuplicatedCandidates(
                rbind(candidate.solutions$candidates,
                      c(2, somatic.purity, NA, 2)))
        }
    }

    if (nrow(candidate.solutions$candidates) > max.candidate.solutions) {
        # test the best solutions and everything close to diploid
        idx.keep <- unique(c(seq_len(max.candidate.solutions), which(
            candidate.solutions$candidates$tumor.ploidy >= max(1.5, min.ploidy) &
            candidate.solutions$candidates$tumor.ploidy <= min(max.ploidy, 2.6) &
            candidate.solutions$candidates$purity >= min(test.purity) &
            candidate.solutions$candidates$purity <= max(test.purity)
            )))
        candidate.solutions$candidates <- candidate.solutions$candidates[idx.keep,
            ]
    }
    flog.info(paste(strwrap(paste("Local optima:\n", paste(round(candidate.solutions$candidates$purity,
        digits = 2), round(candidate.solutions$candidates$ploidy, digits = 2),
        sep = "/", collapse = ", "))), collapse = "\n"))

    simulated.annealing <- TRUE

    .optimizeSolution <- function(cpi) {
        max.attempts <- 4
        attempt <- 0
        while (attempt < max.attempts) {
            attempt <- attempt + 1
            total.ploidy <- candidate.solutions$candidates$ploidy[cpi]
            p <- candidate.solutions$candidates$purity[cpi]
            # optimize purity within +/- 0.2 in the first attempt, then just
            # try to match the ploidy
            idxLocal <- which(abs(test.purity-p) < (0.1+attempt/10))
            test.purity.local <- test.purity[idxLocal]
            prior.purity.local <- prior.purity[idxLocal]

            flog.info("Testing local optimum %i/%i at purity %.2f and total ploidy %.2f...",
                cpi, nrow(candidate.solutions$candidates), p, total.ploidy)

            subclonal <- rep(FALSE, nrow(seg))
            old.llik <- -1
            cnt.llik.equal <- 0
            C.likelihood <- matrix(ncol = length(test.num.copy) + 1, nrow = nrow(seg))
            colnames(C.likelihood) <- c(test.num.copy, "Subclonal")
            sd.seg.orig <- sd.seg
            for (iter in seq_len(iterations)) {
                flog.debug("Iteration %i, Purity %.2f, total ploidy %.2f, ploidy %.2f, likelihood %.3f, offset %.3f error %.3f",
                    iter, p, total.ploidy,  weighted.mean(C, li), llik, mean(log.ratio.offset), sd.seg)
                # test for convergence
                if (abs(old.llik - llik) < 0.0001) {
                  cnt.llik.equal <- cnt.llik.equal + 1
                }
                
                old.llik <- llik
                if (cnt.llik.equal > 3)
                  break
                subclonal.f <- length(unlist(exon.lrs[subclonal])) / length(unlist(exon.lrs))
                # should not happen, but sometimes does for very unlikely local optima.
                if (subclonal.f > 0) flog.debug("Subclonal %.2f", subclonal.f)
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
                  px.rij <- lapply(test.purity.local, function(px) vapply(which(!is.na(C)),
                    function(i) .calcLlikSegment(subclonal = subclonal[i], lr = exon.lrs[[i]] +
                      log.ratio.offset[i], sd.seg = sd.seg, p = px, Ci = C[i], total.ploidy = total.ploidy,
                      max.exon.ratio = max.exon.ratio), double(1)))
                  px.rij.s <- vapply(px.rij, sum, na.rm = TRUE, double(1)) +
                      log(prior.purity.local)

                  if (simulated.annealing) 
                    px.rij.s <- px.rij.s * exp(iter / 4)

                  px.rij.s <- exp(px.rij.s - max(px.rij.s))
                  # Gibbs sample purity
                  p <- test.purity.local[min(which(runif(n = 1, min = 0, max = sum(px.rij.s)) <=
                    cumsum(px.rij.s)))]
                  total.ploidy <- p * (sum(li * (C))) / sum(li) + (1 - p) * 2
                  # Gibbs sample offset
                  if (iter > 2) {
                    log.ratio.offset <- .sampleOffset(subclonal, seg, exon.lrs, sd.seg,
                      p, C, total.ploidy, max.exon.ratio, simulated.annealing, iter,
                      log.ratio.calibration)
                    sd.seg <- .sampleError(subclonal, seg, exon.lrs, sd.seg.orig,
                      p, C, total.ploidy, max.exon.ratio, simulated.annealing, iter,
                      log.ratio.calibration, log.ratio.offset)
                   }
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
                    li[i] * Ci) / sum(li), double(1))

                  # set probability to zero if ploidy is not within requested range
                  log.prior.ploidy <- log(ifelse(ploidy < min.ploidy | ploidy > max.ploidy,
                    0, 1))
                  if (iter > 1)
                    p.rij <- p.rij + log.prior.ploidy

                  liChr <- li
                  idxOtherChrs <- which(seg$chrom != seg$chrom[i])
                  # with small panels, we cannot avoid too long segments and thus
                  # allow chromosomes with lots of deletions in this case
                  if (sum(seg$num.mark[-idxOtherChrs], na.rm = TRUE) > 500) {
                    liChr[idxOtherChrs] <- 0
                  }
                  frac.homozygous.loss <- vapply(test.num.copy, function(Ci) (sum(liChr[-i] *
                    ifelse(C[-i] < 0.5, 1, 0)) + liChr[i] * ifelse(Ci < 0.5, 1, 0)) / sum(liChr),
                    double(1))
                  if (li[i] > max.homozygous.loss[2] && test.num.copy[1] < 1) {
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
                    p.rij <- p.rij * exp(iter / 4)

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
                  opt.C <- (2^(seg$seg.mean + log.ratio.offset) * total.ploidy) / p -
                    ((2 * (1 - p)) / p)
                  opt.C[opt.C < 0] <- 0
                  # if the sub-clonal event is a homozygous loss, we might reach our limit
                  if (id > length(test.num.copy) && opt.C[i] <= 0 &&
                      is.infinite(min(log.prior.homozygous.loss, na.rm = TRUE))) {
                    C[i] <- 1
                    subclonal[i] <- FALSE
                  } else if (id > length(test.num.copy)) {
                    # optimal non-integer copy number
                    C[i] <- opt.C[i]
                    subclonal[i] <- TRUE
                  } else if (test.num.copy[id] == max(test.num.copy) &&
                             opt.C[i] > test.num.copy[id] + 1.5) {
                    # rounded ideal integer copy number is max copy number + 2
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
        list(log.likelihood = llik, purity = p, ploidy = weighted.mean(C, li),
             ML.C = C, opt.C = opt.C, ML.Subclonal = subclonal, total.ploidy = total.ploidy,
             seg = seg.adjusted, C.likelihood = C.likelihood, fraction.subclonal = subclonal.f,
             fraction.homozygous.loss = sum(li[which(C < 0.01)]) / sum(li),
             log.ratio.offset = log.ratio.offset,
             log.ratio.sdev = sd.seg, SA.iterations = iter, candidate.id = cpi)
    }

    .fitSolution <- function(sol) {
        SNV.posterior <- NULL
        p <- sol$purity
        if (!is.null(vcf.file)) {
            flog.info("Fitting variants with %s model for local optimum %i/%i...",
                model, sol$candidate.id, nrow(candidate.solutions$candidates))
        }
        if (sol$fraction.subclonal > max.non.clonal) {
            .stopRuntimeError(".fitSolution received high subclonal solution.")
        }

        gene.calls <- .getGeneCalls(sol$seg, tumor, log.ratio, fun.focal,
            args.focal, chr.hash)
        
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
                .fitSNVp <- function(px, cont.rate = prior.contamination) {
                    flog.info("Fitting variants for purity %.2f, tumor ploidy %.2f and contamination %.2f.",
                        px, weighted.mean(sol$ML.C, sol$seg$size), cont.rate)
                  
                  .calcSNVLLik(vcf, tumor.id.in.vcf, ov, px, test.num.copy,
                    sol$C.likelihood, sol$ML.C, sol$opt.C, median.C = median(rep(sol$ML.C, sol$seg$num.mark)),
                    snv.model = model,
                    snv.lr, sampleid, cont.rate = cont.rate, prior.K = prior.K,
                    max.coverage.vcf = max.coverage.vcf, non.clonal.M = non.clonal.M,
                    model.homozygous = model.homozygous, error = error,
                    max.mapping.bias = max.mapping.bias, max.pon = max.pon,
                    min.variants.segment = min.variants.segment)
                }

                res.snvllik <- lapply(tp[1], .fitSNVp)
                if (post.optimize && length(tp) > 1) {
                    GoF <- .getGoF(list(SNV.posterior = res.snvllik[[1]]))
                    idx <- 1
                    if (GoF < (min.gof - 0.05) && speedup.heuristics > 0) {
                        flog.info("Poor goodness-of-fit (%.3f). Skipping post-optimization.", GoF)
                    } else if (.calcFractionBalanced(res.snvllik[[1]]$posteriors) > 0.8 &&
                        weighted.mean(sol$ML.C, sol$seg$size) > 2.7 && speedup.heuristics > 0) {
                        flog.info("High ploidy solution in highly balanced genome. Skipping post-optimization.")
                    } else if (.isRareKaryotype(weighted.mean(sol$ML.C, sol$seg$size)) && speedup.heuristics > 0) {
                        flog.info("Rare karyotype solution. Skipping post-optimization.")
                    } else {
                        res.snvllik <- c(res.snvllik, lapply(tp[-1], .fitSNVp))
                      px.rij <- lapply(tp, function(px) vapply(which(!is.na(C)), function(i)
                        .calcLlikSegment(subclonal = sol$ML.Subclonal[i],
                        lr = exon.lrs[[i]] + sol$log.ratio.offset[i], sd.seg = sol$log.ratio.sdev, p = px, 
                        Ci = sol$ML.C[i], total.ploidy = px * (sum(sol$seg$size * sol$ML.C))/sum(sol$seg$size) + (1 - 
                          px) * 2, max.exon.ratio = max.exon.ratio), double(1)))

                      px.rij.s <- vapply(px.rij, sum, na.rm = TRUE, double(1)) + 
                          log(pp) + vapply(res.snvllik, 
                        function(x) x$llik, double(1))
                      idx <- which.max(px.rij.s)
                  }
                } else {
                  idx <- 1
                }
                p <- tp[idx]
                flog.info("Optimized purity: %.2f", p)
                SNV.posterior <- res.snvllik[[idx]]
                list(p = p, SNV.posterior = SNV.posterior)
            }
            fitRes <- .fitSNV(tp, pp)
            sol$purity <- fitRes$p
            SNV.posterior <- fitRes$SNV.posterior
        }

        c(sol, list(
            C.posterior = data.frame(sol$C.likelihood/rowSums(sol$C.likelihood), ML.C =
            sol$ML.C, Opt.C = sol$opt.C, ML.Subclonal = sol$ML.Subclonal),
            SNV.posterior = SNV.posterior,
             fraction.balanced = 
            .calcFractionBalanced(SNV.posterior$posteriors),
            gene.calls = gene.calls,
            failed = FALSE))
    }

    # assign integer copy numbers to segments using SA 
    if (is.null(BPPARAM)) {
        results <- lapply(seq_len(nrow(candidate.solutions$candidates)), .optimizeSolution)
    } else {
        results <- BiocParallel::bplapply(seq_len(nrow(candidate.solutions$candidates)), 
                                          .optimizeSolution, BPPARAM = BPPARAM)
    }
    # rank and then delete lower ranked similar solution
    results <- .rankResults(results)
    nBefore <- length(results)
    results <- .filterDuplicatedResults(results, purity.cutoff = 0.05)
    # bring back in original order for progress output
    results <- results[order(vapply(results, function(sol) sol$candidate.id, integer(1)))]

    if (length(results) < nBefore) { 
        flog.info("Skipping %i solutions that converged to the same optima.",
                  nBefore - length(results))
    }
    idxFailed <- vapply(results, function(sol)
                        sol$fraction.subclonal > max.non.clonal, logical(1))
    if (sum(is.na(idxFailed))) .stopRuntimeError("NAs in fraction.subclonal.")
    if (sum(idxFailed)) {
        tmp <- sapply(which(idxFailed), function(i) paste0(results[[i]]$purity, "/", round(results[[i]]$ploidy, digits = 2)))
        flog.info("Skipping %i solutions exceeding max.non.clonal (%.2f): %s (purity/tumor ploidy)",
                  sum(idxFailed), max.non.clonal, paste(tmp, collapse=", "))
    }
    results <- results[which(!idxFailed)]

    # fit SNVs
    if (!is.null(vcf.file) && .countVariants(vcf) < min.variants) {
        flog.warn("Insufficient number of variants %i passing filters. There is an issue with your data.",
                  .countVariants(vcf))
    } else {
        if (is.null(BPPARAM) || is.null(vcf.file)) {
            results <- lapply(results, .fitSolution)
        } else {
            results <- BiocParallel::bptry(BiocParallel::bplapply(results,
                                              .fitSolution, BPPARAM = BPPARAM))

            if (!all(BiocParallel::bpok(results))) {
                .stopRuntimeError("Error in parallel variant fitting: ",
                    tail(attr(results[[which(!BiocParallel::bpok(results))]], "traceback")),
                    "\nBefore reporting an issue, try again with fewer cores and/or more memory."
                )
            }
        }
    }
    results <- .rankResults(results)
    results <- .filterDuplicatedResults(results)

    # If necessary, try to estimate contamination
    #--------------------------------------------------------------------------
    cont.rate <- prior.contamination

    for (i in seq_along(results)) {
        if (grepl("CONTAMINATION", vcf.filtering$flag_comment)) {
            cont.rate <- .estimateContamination(
                results[[i]]$SNV.posterior$posteriors,
                max.mapping.bias)
            if (cont.rate > prior.contamination) {
                flog.info("Initial guess of contamination rate: %.3f", cont.rate)
            }
        }
        ## optimize contamination. we just re-run the fitting if necessary (
        if (grepl("CONTAMINATION", vcf.filtering$flag_comment)) {
            if (cont.rate > prior.contamination) {
                flog.info("Optimizing contamination rate of optimum %i/%i...",
                    i, length(results))

                res.snvllik <-
                    .calcSNVLLik(vcf, tumor.id.in.vcf,
                          ov, results[[i]]$purity, test.num.copy, results[[i]]$C.likelihood,
                          results[[i]]$C.posterior$ML.C,
                          results[[i]]$C.posterior$Opt.C,
                          median.C = median(rep(results[[i]]$seg$C, results[[i]]$seg$num.mark)),
                          snv.model = model,
                          snv.lr, sampleid, cont.rate = cont.rate, prior.K = prior.K,
                          max.coverage.vcf = max.coverage.vcf, non.clonal.M = non.clonal.M,
                          model.homozygous = model.homozygous, error = error,
                          max.mapping.bias = max.mapping.bias, max.pon = max.pon,
                          min.variants.segment = min.variants.segment)
                results[[i]]$SNV.posterior <- res.snvllik
            }
            cont.rate <- .estimateContamination(
                        results[[i]]$SNV.posterior$posteriors,
                                    max.mapping.bias)
                    flog.info("Optimized contamination rate: %.3f", cont.rate)
            results[[i]]$SNV.posterior$posterior.contamination <- cont.rate
        }
    }
    # re-order after contamination optimization
    if (length(results) > 1) results <- .rankResults(results)

    if (length(results) &&
        !is.null(results[[1]]$SNV.posterior$posterior.contamination) &&
        results[[1]]$SNV.posterior$posterior.contamination < 0.001 &&
        !is.na(vcf.filtering$flag_comment) &&
        vcf.filtering$flag_comment == "POTENTIAL SAMPLE CONTAMINATION") {
        vcf.filtering$flag <- FALSE
        vcf.filtering$flag_comment <- ""
    }

    results <- .flagResults(results, max.non.clonal = max.non.clonal, max.logr.sdev = max.logr.sdev,
        logr.sdev = sd.seg, max.segments = max.segments, min.gof = min.gof, flag = vcf.filtering$flag,
        flag_comment = vcf.filtering$flag_comment, dropout = dropoutWarning,
        use.somatic.status = args.filterVcf$use.somatic.status,
        model.homozygous = model.homozygous)

    if (length(results) < 1) {
        flog.warn("Could not find valid purity and ploidy solution.")
    }

    if (grepl("NON-ABERRANT", results[[1]]$flag_comment)) {
        flog.warn("No copy number alterations found, purity estimate unreliable.")
        results <- .findClosestSolution(results,
            purity = min(test.purity), ploidy = 2, ploidy.div = 1)
    }
    .logFooter()
    rho <- if (is.null(vcf)) NULL else info(vcf)[[paste0(vcf.field.prefix, "MBRHO")]]
    list(
        candidates = candidate.solutions,
        results = results,
        exon.lrs = exon.lrs,
        input = list(tumor = tumor.coverage.file, normal = normal.coverage.file,
            log.ratio = GRanges(normal[, 1], on.target = normal$on.target, log.ratio = log.ratio),
            log.ratio.sdev = sd.seg, mapd = mapd, vcf = vcf, sampleid = sampleid,
            test.num.copy = test.num.copy,
            sex = sex, sex.vcf = sex.vcf, chr.hash = chr.hash, centromeres = centromeres,
            mapping.bias.rho = if (is.null(rho) || all(is.na(rho))) NULL else mean(rho, na.rm = TRUE),
            vcf.field.prefix = vcf.field.prefix,
            interval.file = interval.file,
            args = list(
                filterVcf = args.filterVcf[vapply(args.filterVcf, object.size, double(1)) < 1000],
                filterIntervals = args.filterIntervals[vapply(args.filterIntervals, object.size, double(1)) < 1000])
        ),
        version = package.version("PureCN")
    )
}
