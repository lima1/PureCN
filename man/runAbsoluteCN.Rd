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
\usage{runAbsoluteCN(gatk.normal.file = NULL, gatk.tumor.file, 
    log.ratio = NULL, seg.file = NULL, seg.file.sdev = 0.4, 
    vcf.file = NULL, genome, centromeres = NULL, sex = c("?", 
        "F", "M", "diploid"), fun.filterVcf = filterVcfMuTect, 
    args.filterVcf = list(), fun.setPriorVcf = setPriorVcf, 
    args.setPriorVcf = list(), fun.segmentation = segmentationCBS, 
    args.segmentation = list(), fun.focal = findFocal, 
    args.focal = list(), sampleid = NULL, min.ploidy = 1, 
    max.ploidy = 6, test.num.copy = 0:7, test.purity = seq(0.15, 
        0.95, by = 0.01), prior.purity = rep(1, length(test.purity))/length(test.purity), 
    max.candidate.solutions = 20, candidates = NULL, 
    coverage.cutoff = 15, max.non.clonal = 0.2, max.homozygous.loss = 0.1, 
    iterations = 30, log.ratio.calibration = 0.25, 
    remove.off.target.snvs = TRUE, gc.gene.file = NULL, 
    filter.lowhigh.gc.exons = 0.001, max.dropout = c(0.95, 
        1.1), filter.targeted.base = 4, max.logr.sdev = 0.75, 
    max.segments = 200, prior.K = 0.999, chr.hash = NULL, 
    plot.cnv = TRUE, verbose = TRUE, post.optimize = FALSE, 
    ...)}
\arguments{
  \item{gatk.normal.file}{GATK coverage file of normal control (optional if 
log.ratio is provided - then it will be only used to filter low coverage 
exons). 
Should be already GC-normalized with \code{\link{correctCoverageBias}}. 
Needs to be either a file name or data read with the 
\code{\link{readCoverageGatk}} function.}
  \item{gatk.tumor.file}{GATK coverage file of tumor. 
Should be already GC-normalized with \code{\link{correctCoverageBias}}. 
Needs to be either a file name or data read with the 
\code{\link{readCoverageGatk}} function.}
  \item{log.ratio}{Copy number log-ratios for all exons in the coverage files. 
If NULL, calculated based on coverage files.}
  \item{seg.file}{Segmented data. Optional, to support matched SNP6 data. 
If NULL, use coverage files or log.ratio to segment the data.  }
  \item{seg.file.sdev}{If seg.file provided, the log-ratio standard deviation, 
used to model likelihood of sub-clonal copy number events.}
  \item{vcf.file}{VCF file, tested with MuTect output files.  Optional, but
typically needed to select between local optima of similar likelihood. Can 
also be a CollapsedVCF, read with the \code{readVcf} function. 
Requires a DB info flag for dbSNP membership. The default 
\code{fun.setPriorVcf} function will also
look for a Cosmic.CNT slot, containing the hits in the COSMIC database.
Again, do not expect very useful results without a VCF file.}
  \item{genome}{Genome version, required for the \code{readVcf} function.}
  \item{centromeres}{Data.frame with centromere positions in first three columns. 
If NULL, use pre-stored positions for genome versions hg18, hg19 and hg38. }
  \item{sex}{Sex of sample. If ?, detect using \code{\link{getSexFromCoverage}} 
function and default parameters. 
Default parameters might not work well with every assay and 
might need to be tuned. If set to diploid, then PureCN will assume all 
chromosomes are diploid and will not try to detect sex. }
  \item{fun.filterVcf}{Function for filtering variants. Expected output is a 
list with elements vcf (CollapsedVCF), flag (TRUE/FALSE) and flag_comment 
(string). The flags will be added to the output data and can be used to 
warn users, for example when samples look too noisy. Default filter will 
remove variants flagged by MuTect, but will keep germline variants. If 
ran in matched normal mode, it will by default use somatic status of 
variants and filter non-somatic calls with allelic fraction significantly 
different from 0.5 in normal. }
  \item{args.filterVcf}{Arguments for variant filtering function. Arguments 
vcf, tumor.id.in.vcf, coverage.cutoff and verbose are required in the 
filter function and are automatically set (do NOT set them here again).}
  \item{fun.setPriorVcf}{Function to set prior for somatic status for each 
variant in the VCF.}
  \item{args.setPriorVcf}{Arguments for somatic prior function.}
  \item{fun.segmentation}{Function for segmenting the copy number log-ratios. 
Expected return value is a list with elements seg (the segmentation) and 
size (the size in bp for all segments).}
  \item{args.segmentation}{Arguments for segmentation function. Arguments 
normal, tumor, log.ratio, plot.cnv, coverage.cutoff, sampleid, vcf, 
tumor.id.in.vcf, verbose are required in the segmentation function and 
automatically set (do NOT set them here again).}
  \item{fun.focal}{Function for identifying focal amplifications.}
  \item{args.focal}{Arguments for focal amplification function.}
  \item{sampleid}{Sample id, provided in output files etc.}
  \item{min.ploidy}{Minimum ploidy to be considered.}
  \item{max.ploidy}{Maximum ploidy to be considered.}
  \item{test.num.copy}{Copy numbers tested in the grid search. Note that focal 
amplifications can have much higher copy numbers, but they will be labeled
as subclonal (because they do not fit the integer copy numbers).}
  \item{test.purity}{Considered tumor purity values. }
  \item{prior.purity}{Priors for purity if they are available. Only change 
when you know what you are doing.}
  \item{max.candidate.solutions}{Number of local optima considered in optimization 
and variant fitting steps. If there are too many local optima, it will use
specified number of top candidate solutions, but will also include all 
optima close to diploid, because silent genomes have often lots of local 
optima.}
  \item{candidates}{Candidates to optimize from a previous run 
(\code{return.object$candidates}). 
If NULL, do 2D grid search and find local optima. }
  \item{coverage.cutoff}{Minimum exon coverage in both normal and tumor. Exons
with lower coverage are ingored. The cutoff choice depends on the expected
purity and overall coverage. High purity samples might need a lower cutoff
to call homozygous deletions. If an \code{exon.weight.file} (below) is NOT 
specified, it is recommended to set a higher cutoff (e.g. 20) to remove 
noise from unreliable exon measurements. }
  \item{max.non.clonal}{Maximum genomic fraction assigned to a subclonal copy 
number state.}
  \item{max.homozygous.loss}{Maximum genomic fraction assigned to homozygous loss. 
This is set to a fairly high default value to not exclude correct
solutions, especially in noisy segmentations. }
  \item{iterations}{Maximum number of iterations in the Simulated Annealing copy 
number fit optimization.}
  \item{log.ratio.calibration}{Re-calibrate log-ratios in the window 
sd(log.ratio)*log.ratio.calibration.}
  \item{remove.off.target.snvs}{If set to a true value, will remove all SNVs outside the 
covered regions.}
  \item{gc.gene.file}{A mapping file that assigns GC content and gene symbols 
to each exon in the coverage files. Used for generating gene-level calls. 
First column in format CHR:START-END. Second column GC content (0 to 1). 
Third column gene symbol. This file
can be generated with the \sQuote{GATK GCContentByInterval} tool or 
with the \code{\link{calculateGCContentByInterval}} function.}
  \item{filter.lowhigh.gc.exons}{Quantile q (defines lower q and upper 1-q) 
for removing exons with outlier GC profile. Assuming that GC correction 
might not have been worked on those. Requires gc.gene.file.}
  \item{max.dropout}{Measures GC bias as ratio of coverage in AT-rich (GC < 0.5) 
versus GC-rich regions (GC >= 0.5). High drop-out might indicate that
data was not GC-normalized or that the sample quality might be 
insufficient. Requires gc.gene.file.}
  \item{filter.targeted.base}{Exclude exons with targeted base (size) smaller 
than this cutoff. This is useful when the same interval file was used to
calculate GC content. For such small exons, the GC content is likely 
very different from the true GC content of the probes.}
  \item{max.logr.sdev}{Flag noisy samples with segment log-ratio standard deviation 
larger than this. Assay specific and needs to be calibrated.}
  \item{max.segments}{Flag noisy samples with a large number of segments. Assay 
specific and needs to be calibrated.}
  \item{prior.K}{This defines the prior probability that the multiplicity of
a SNV corresponds to either the maternal or the paternal copy
number (for somatic variants additionally a multiplicity of 1).
For perfect segmentations, this value would be 1; values smaller
than 1 thus may provide some robustness against segmentation errors.}
  \item{chr.hash}{Mapping of non-numerical chromsome names to numerical names
(e.g. chr1 to 1, chr2 to 2, etc.). If NULL, assume chromsomes
are properly ordered.   }
  \item{plot.cnv}{Generate segmentation plots.}
  \item{verbose}{Verbose output.}
  \item{post.optimize}{Optimize purity using final SCNA-fit and SNVs. This 
might take a long time when lots of SNVs need to be fitted, but will 
typically result in a slightly more accurate purity, especially for rather
silent genomes or very low purities. Otherwise, it will just use the 
purity determined via the SCNA-fit.}
  \item{\dots}{Additional parameters passed to the segmentation function.}
}

\value{A list with elements
\item{candidates}{Results of the grid search.}
\item{results}{All local optima, sorted by final rank.}
\item{input}{The input data.}}

\author{Markus Riester}



\seealso{\code{\link{correctCoverageBias} \link{segmentationCBS}}}
\examples{
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
}
