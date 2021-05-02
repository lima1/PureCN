#' GATK4 ModelSegments segmentation function
#' 
#' A wrapper for GATK4s ModelSegmentation function, useful when normalization
#' is performed with other tools than GATK4, for example PureCN.
#' This function is called via the
#' \code{fun.segmentation} argument of \code{\link{runAbsoluteCN}}.  The
#' arguments are passed via \code{args.segmentation}.
#' 
#' 
#' @param normal Coverage data for normal sample. Ignored in this function.
#' @param tumor Coverage data for tumor sample.
#' @param log.ratio Copy number log-ratios, one for each exon in coverage file.
#' @param seg If segmentation was provided by the user, this data structure
#' will contain this segmentation. Useful for minimal segmentation functions.
#' Otherwise PureCN will re-segment the data. This segmentation function
#' ignores this user provided segmentation.
#' @param vcf Optional \code{CollapsedVCF} object with germline allelic ratios.
#' @param tumor.id.in.vcf Id of tumor in case multiple samples are stored in
#' VCF.
#' @param normal.id.in.vcf Id of normal in in VCF. Currently not used.
#' @param prune.hclust.h Ignored in this function.
#' @param prune.hclust.method Ignored in this function.
#' @param additional.cmd.args \code{character(1)}. By default,
#' \code{ModelSegments} is called with default parameters. Provide additional 
#' arguments here.
#' @param chr.hash Not needed here since \code{ModelSegments} does not
#' require numbered chromosome names.
#' @param ... Currently unused arguments provided to other segmentation
#' functions.
#' @return \code{data.frame} containing the segmentation.
#' @author Markus Riester
#'
#' @seealso \code{\link{runAbsoluteCN}}
#' @examples
#' 
#' normal.coverage.file <- system.file("extdata", "example_normal_tiny.txt", 
#'     package="PureCN")
#' tumor.coverage.file <- system.file("extdata", "example_tumor_tiny.txt", 
#'     package="PureCN")
#' vcf.file <- system.file("extdata", "example.vcf.gz",
#'     package="PureCN")
#' 
#' # The max.candidate.solutions, max.ploidy and test.purity parameters are set to
#' # non-default values to speed-up this example.  This is not a good idea for real
#' # samples.
#'\dontrun{
#'  ret <-runAbsoluteCN(normal.coverage.file=normal.coverage.file, 
#'      tumor.coverage.file=tumor.coverage.file, vcf.file=vcf.file, 
#'      sampleid="Sample1",  genome="hg19",
#'      fun.segmentation = segmentationGATK4, max.ploidy=4,
#'      args.segmentation = list(additional.cmd.args = "--gcs-max-retries 19"),
#'      test.purity=seq(0.3,0.7,by=0.05), max.candidate.solutions=1)
#'}
#' 
#' @importFrom utils compareVersion
#' @export segmentationGATK4
segmentationGATK4 <- function(normal, tumor, log.ratio, seg, 
    vcf = NULL, tumor.id.in.vcf = 1, normal.id.in.vcf = NULL,
    prune.hclust.h = NULL, prune.hclust.method = NULL,
    additional.cmd.args = "",
    chr.hash = NULL, ...) {

    min.version <- "4.1.7.0"
    if (.checkGATK4Version(min.version) < 0) {
        .stopUserError("Unable to find a recent (>= %s) gatk binary.",
            min.version)
    }    
    output.vcf.file <- NULL    
    if (!is.null(vcf)) {
        output.vcf.file <- tempfile(fileext = ".tsv")
        .writeAllelicCountsFileGatk(vcf, tumor.id.in.vcf, output.vcf.file)
    } 
    if (length(log.ratio) != length(tumor)) {
        .stopRuntimeError("tumor and log.ratio do not align in segmentationGATK4.")
    }
    output.lr.file <- tempfile(fileext = ".tsv")
    # following .writeLogRatioFileGATK4 needs the GRanges information, so
    # simply update it (should be the same though!)
    tumor$log.ratio <- log.ratio
    .writeLogRatioFileGATK4(list(vcf=vcf, log.ratio = tumor), tumor.id.in.vcf, output.lr.file)
    output.dir <- tempdir()
    sampleid <- .getSampleIdFromVcf(vcf, tumor.id.in.vcf)
    args <- paste("ModelSegments",
              "--allelic-counts", output.vcf.file,
              "--denoised-copy-ratios", output.lr.file,
              "--output-prefix", "tumor",
              additional.cmd.args,
              "-O", output.dir)
    output <- try(system2("gatk", args, stderr = TRUE, stdout = TRUE))
    if (is(output, "try-error")) {
        .stopUserError("GATK4 ModelSegments runtime error: ",
            output, "\nArguments: ", args)
    }
    seg <- readSegmentationFile(file.path(output.dir, "tumor.modelFinal.seg"), sampleid)
    idx.enough.markers <- seg$num.mark > 1
    rownames(seg) <- NULL
    file.remove(output.vcf.file)
    file.remove(output.lr.file)
    attr(seg, "ModelSegments.Output") <- output
    attr(seg, "ModelSegments.Args") <- args
    seg[idx.enough.markers,]
}

.getSampleIdFromVcf <- function(vcf, tumor.id.in.vcf) {
    if (is(tumor.id.in.vcf, "character") &&
        tumor.id.in.vcf %in% samples(header(vcf))) return(tumor.id.in.vcf)
    if (is(tumor.id.in.vcf, "numeric") &&
        length(samples(header(vcf))) <= tumor.id.in.vcf) {
        return(samples(header(vcf))[tumor.id.in.vcf])
    }
    return(NA)    
}    

.checkGATK4Version <- function(min.version) {
    output <- try(system2("gatk", "--version", stdout=TRUE), silent = TRUE)
    if (is(output, "try-error")) {
        flog.warn("Cannot find gatk binary in path.")
        return(-1)
    }
    prefix <- "The Genome Analysis Toolkit (GATK) v"    
    compareVersion(gsub(prefix, "", output[grep(prefix, output, fixed = TRUE)], fixed = TRUE), min.version)
}    
