\name{runAbsoluteCN}
\alias{runAbsoluteCN}
\title{Run PureCN implementation of ABSOLUTE}
\description{This function takes as input tumor and normal control coverage and 
allelic fractions of germline variants and somatic mutations.
Coverage data is provided in GATK DepthOfCoverage format, allelic fraction
in VCF format (e.g. obtained by MuTect). Normal control does not need to 
be matched (from the same patient). In case VCF does not contain somatic 
status, it should contain dbSNP and optionally COSMIC annotation.
Returns purity and ploidy combinations, sorted by likelihood score.
Provides copy number and LOH data, by both gene and genomic region. }
\usage{runAbsoluteCN(normal.coverage.file = NULL, tumor.coverage.file = NULL, 
    log.ratio = NULL, seg.file = NULL, seg.file.sdev = 0.4, 
    vcf.file = NULL, genome, centromeres = NULL, sex = c("?", 
        "F", "M", "diploid"), fun.filterVcf = filterVcfMuTect, 
    args.filterVcf = list(), fun.setPriorVcf = setPriorVcf, 
    args.setPriorVcf = list(), fun.setMappingBiasVcf = setMappingBiasVcf, 
    args.setMappingBiasVcf = list(), fun.filterTargets = filterTargets, 
    args.filterTargets = list(), fun.segmentation = segmentationCBS, 
    args.segmentation = list(), fun.focal = findFocal, 
    args.focal = list(), sampleid = NULL, min.ploidy = 1, 
    max.ploidy = 6, test.num.copy = 0:7, test.purity = seq(0.15, 
        0.95, by = 0.01), prior.purity = NULL, prior.K = 0.999, 
    prior.contamination = 0.01, max.candidate.solutions = 20, 
    candidates = NULL, coverage.cutoff = NULL, min.coverage = 15, 
    max.coverage.vcf = 300, max.non.clonal = 0.2, max.homozygous.loss = 0.1, 
    iterations = 30, log.ratio.calibration = 0.25, 
    remove.off.target.snvs = TRUE, gc.gene.file = NULL, 
    max.dropout = c(0.95, 1.1), max.logr.sdev = 0.75, 
    max.segments = 200, plot.cnv = TRUE, cosmic.vcf.file = NULL, 
    post.optimize = FALSE, verbose = TRUE, gatk.tumor.file = NULL, 
    gatk.normal.file = NULL, ...)}
\arguments{
  \item{normal.coverage.file}{GATK coverage file of normal control (optional if 
log.ratio is provided - then it will be only used to filter low coverage 
exons). 
Should be already GC-normalized with \code{\link{correctCoverageBias}}. 
Needs to be either a file name or data read with the 
\code{\link{readCoverageGatk}} function.}
  \item{tumor.coverage.file}{GATK coverage file of tumor. If \code{NULL}, requires \code{seg.file}
and an interval file via \code{gc.gene.file}.
Should be already GC-normalized with \code{\link{correctCoverageBias}}. 
Needs to be either a file name or data read with the 
\code{\link{readCoverageGatk}} function.}
  \item{log.ratio}{Copy number log-ratios for all exons in the coverage files. 
If \code{NULL}, calculated based on coverage files.}
  \item{seg.file}{Segmented data. Optional, to support matched SNP6 data. 
If \code{NULL}, use coverage files or \code{log.ratio} to segment the data.}
  \item{seg.file.sdev}{If \code{seg.file} provided, the log-ratio standard deviation, 
used to model likelihood of sub-clonal copy number events.}
  \item{vcf.file}{VCF file, tested with \sQuote{MuTect 1} output files.  Optional, but
typically needed to select between local optima of similar likelihood. Can 
also be a \code{CollapsedVCF}, read with the \code{readVcf} function. 
Requires a DB info flag for dbSNP membership. The default 
\code{fun.setPriorVcf} function will also
look for a Cosmic.CNT slot (see \code{cosmic.vcf.file}), containing 
the hits in the COSMIC database.
Again, do not expect very useful results without a VCF file.}
  \item{genome}{Genome version, for example hg19.}
  \item{centromeres}{A \code{data.frame} with centromere positions in first three columns. 
If \code{NULL}, use pre-stored positions for genome versions 
hg18, hg19 and hg38. }
  \item{sex}{Sex of sample. If \code{?}, detect using \code{\link{getSexFromCoverage}} 
function and default parameters. 
Default parameters might not work well with every assay and 
might need to be tuned. If set to diploid, then PureCN will assume all 
chromosomes are diploid and will not try to detect sex. }
  \item{fun.filterVcf}{Function for filtering variants. Expected output is a 
list with elements \code{vcf} (\code{CollapsedVCF}), 
flag (\code{logical(1)}) and \code{flag_comment} (\code{character(1)}).
The flags will be added to the output data and can be used to 
warn users, for example when samples look too noisy. Default filter will 
remove variants flagged by MuTect, but will keep germline variants. If 
ran in matched normal mode, it will by default use somatic status of 
variants and filter non-somatic calls with allelic fraction significantly 
different from 0.5 in normal. Defaults to \code{\link{filterVcfMuTect}}, 
which in turn also calls \code{\link{filterVcfBasic}}.}
  \item{args.filterVcf}{Arguments for variant filtering function. Arguments 
\code{vcf}, \code{tumor.id.in.vcf}, \code{min.coverage} and 
\code{verbose} are required in the 
filter function and are automatically set (do NOT set them here again).}
  \item{fun.setPriorVcf}{Function to set prior for somatic status for each 
variant in the VCF. Defaults to \code{\link{setPriorVcf}}.}
  \item{args.setPriorVcf}{Arguments for somatic prior function.}
  \item{fun.setMappingBiasVcf}{Function to set mapping bias for each 
variant in the VCF. Defaults to \code{\link{setMappingBiasVcf}}.}
  \item{args.setMappingBiasVcf}{Arguments for mapping bias function.}
  \item{fun.filterTargets}{Function for filtering low-quality targets in the coverage 
files. Needs to return a \code{logical} vector whether an interval
should be used for segmentation. Defaults to \code{\link{filterTargets}}.}
  \item{args.filterTargets}{Arguments for target filtering function. Arguments
\code{log.ratio}, \code{tumor}, \code{gc.data}, \code{seg.file},
\code{verbose} are required and automatically set (do NOT set 
them here again).}
  \item{fun.segmentation}{Function for segmenting the copy number log-ratios. 
Expected return value is a \code{data.frame} representation of the 
segmentation. Defaults to \code{\link{segmentationCBS}}.}
  \item{args.segmentation}{Arguments for segmentation function. Arguments 
\code{normal}, \code{tumor}, \code{log.ratio}, \code{plot.cnv}, 
\code{min.coverage}, \code{sampleid}, \code{vcf}, 
\code{tumor.id.in.vcf}, \code{verbose} are required in the 
segmentation function and automatically set (do NOT set them here again).}
  \item{fun.focal}{Function for identifying focal amplifications. Defaults to 
\code{\link{findFocal}}.}
  \item{args.focal}{Arguments for focal amplification function.}
  \item{sampleid}{Sample id, provided in output files etc.}
  \item{min.ploidy}{Minimum ploidy to be considered.}
  \item{max.ploidy}{Maximum ploidy to be considered.}
  \item{test.num.copy}{Copy numbers tested in the grid search. Note that focal 
amplifications can have much higher copy numbers, but they will be labeled
as subclonal (because they do not fit the integer copy numbers).}
  \item{test.purity}{Considered tumor purity values. }
  \item{prior.purity}{\code{numeric(length(test.purity))} with priors for tested 
purity values. If \code{NULL}, use flat priors. }
  \item{prior.K}{This defines the prior probability that the multiplicity of
a SNV corresponds to either the maternal or the paternal copy
number (for somatic variants additionally to a multiplicity of 1).
For perfect segmentations, this value would be 1; values smaller
than 1 thus may provide some robustness against segmentation errors.}
  \item{prior.contamination}{The prior probability that a known SNP is from a different
individual. }
  \item{max.candidate.solutions}{Number of local optima considered in optimization 
and variant fitting steps. If there are too many local optima, it will use
specified number of top candidate solutions, but will also include all 
optima close to diploid, because silent genomes have often lots of local 
optima.}
  \item{candidates}{Candidates to optimize from a previous run 
(\code{return.object$candidates}). 
If \code{NULL}, do 2D grid search and find local optima. }
  \item{coverage.cutoff}{Deprecated, use \code{min.coverage} instead. }
  \item{min.coverage}{Minimum coverage in both normal and tumor. Targets
with lower coverage are ingored. }
  \item{max.coverage.vcf}{This will set the maximum number of reads in the SNV fitting. 
This is to avoid that small non-reference biases that come apparent only 
at high coverages have a dramatic influence on likelihood scores.}
  \item{max.non.clonal}{Maximum genomic fraction assigned to a subclonal copy 
number state.}
  \item{max.homozygous.loss}{Maximum genomic fraction assigned to homozygous loss. 
This is set to a fairly high default value to not exclude correct
solutions, especially in noisy segmentations. }
  \item{iterations}{Maximum number of iterations in the Simulated Annealing copy 
number fit optimization. Note that this an integer optimization problem
that should converge quickly. Allowed range is 10 to 250.}
  \item{log.ratio.calibration}{Re-calibrate log-ratios in the window 
\code{sd(log.ratio)*log.ratio.calibration}.}
  \item{remove.off.target.snvs}{If set to a true value, will remove all SNVs outside the 
covered regions.}
  \item{gc.gene.file}{A mapping file that assigns GC content and gene symbols 
to each exon in the coverage files. Used for generating gene-level calls. 
First column in format CHR:START-END. Second column GC content (0 to 1). 
Third column gene symbol. This file
can be generated with the \sQuote{GATK GCContentByInterval} tool or 
with the \code{\link{calculateGCContentByInterval}} function.}
  \item{max.dropout}{Measures GC bias as ratio of coverage in AT-rich (GC < 0.5) 
versus GC-rich regions (GC >= 0.5). High drop-out might indicate that
data was not GC-normalized or that the sample quality might be 
insufficient. Requires \code{gc.gene.file}.}
  \item{max.logr.sdev}{Flag noisy samples with segment log-ratio standard deviation 
larger than this. Assay specific and needs to be calibrated.}
  \item{max.segments}{Flag noisy samples with a large number of segments. Assay 
specific and needs to be calibrated.}
  \item{plot.cnv}{Generate segmentation plots.}
  \item{cosmic.vcf.file}{Add a \code{Cosmic.CNT} info field to the provided 
\code{vcf.file} using a VCF file containing the COSMIC database.
The default \code{fun.setPriorVcf} function will give SNVs
found in the COSMIC database a higher prior probability of being somatic.
Not used in likelhood model when matched normal is available in 
\code{vcf.file}. Should be compressed and indexed with bgzip and tabix,
respectively.}
  \item{post.optimize}{Optimize purity using final SCNA-fit and SNVs. This 
might take a long time when lots of SNVs need to be fitted, but will 
typically result in a slightly more accurate purity, especially for rather
silent genomes or very low purities. Otherwise, it will just use the 
purity determined via the SCNA-fit.}
  \item{verbose}{Verbose output.}
  \item{gatk.tumor.file}{Deprecated, use tumor.coverage.file instead.}
  \item{gatk.normal.file}{Deprecated, use normal.coverage.file instead.}
  \item{\dots}{Additional parameters passed to the segmentation function.}
}

\value{A list with elements
\item{candidates}{Results of the grid search.}
\item{results}{All local optima, sorted by final rank.}
\item{input}{The input data.}}

\author{Markus Riester}



\seealso{\code{\link{correctCoverageBias} \link{segmentationCBS}}}
\examples{
normal.coverage.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")
vcf.file <- system.file("extdata", "example_vcf.vcf", 
    package="PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
    package="PureCN")

# The max.candidate.solutions, max.ploidy and test.purity parameters are set to
# non-default values to speed-up this example.  This is not a good idea for real
# samples.
ret <-runAbsoluteCN(normal.coverage.file=normal.coverage.file, 
    tumor.coverage.file=tumor.coverage.file, genome="hg19", vcf.file=vcf.file,
    sampleid='Sample1', gc.gene.file=gc.gene.file,
    max.ploidy=4, test.purity=seq(0.3,0.7,by=0.05), max.candidate.solutions=1)


# If a high-quality segmentation was obtained with third-party tools:
seg.file <- system.file("extdata", "example_seg.txt", 
    package = "PureCN")

# By default, PureCN will re-segment the data, for example to identify
# regions of copy number neutral LOH. If this is not wanted, we can provide
# a minimal segmentation function which just returns the provided one:
funSeg <- function(seg, ...) return(seg)

res <- runAbsoluteCN(seg.file=seg.file, fun.segmentation=funSeg, max.ploidy = 4,
    test.purity = seq(0.3, 0.7, by = 0.05), max.candidate.solutions=1,
    genome='hg19', gc.gene.file=gc.gene.file)
}
