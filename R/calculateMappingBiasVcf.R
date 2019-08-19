#' Calculate Mapping Bias 
#'
#' Function calculate mapping bias for each variant in the provided
#' panel of normals VCF.
#'
#'
#' @param normal.panel.vcf.file Combined VCF file of a panel of normals,
#' reference and alt counts as AD genotype field. Should be compressed and
#' indexed with bgzip and tabix, respectively.
#' @param min.normals Minimum number of normals with heterozygous SNP for
#' calculating position-specific mapping bias. 
#' @param min.normals.betafit Minimum number of normals with heterozygous SNP
#' fitting a beta distribution
#' @param min.median.coverage.betafit Minimum median coverage of normals with
#' heterozygous SNP for fitting a beta distribution
#' @param betafit.coverage.outliers Exclude samples with coverages below or above
#' the specified cutoffs (fractions of coverages median).
#' @param yieldSize See \code{TabixFile}
#' @param genome See \code{readVcf}
#' @return A \code{GRanges} object with mapping bias and number of normal
#' samples with this variant.
#' @author Markus Riester
#' @examples
#'
#' normal.panel.vcf <- system.file("extdata", "normalpanel.vcf.gz", package="PureCN")
#' bias <- calculateMappingBiasVcf(normal.panel.vcf, genome = "h19")
#' saveRDS(bias, "mapping_bias.rds")
#'
#' @importFrom GenomicRanges GRangesList
#' @importFrom fitdistrplus fitdist
#' @export calculateMappingBiasVcf
calculateMappingBiasVcf <- function(normal.panel.vcf.file, min.normals = 2,
                                    min.normals.betafit = 7,
                                    min.median.coverage.betafit = 5,
                                    betafit.coverage.outliers = c(0.25, 4),
                                    yieldSize = 5000, genome) {
    tab <- TabixFile(normal.panel.vcf.file, yieldSize = yieldSize)
    open(tab)
    param <- ScanVcfParam(geno = c("AD"), fixed = "ALT", info = NA)
    cntVar <- 0
    cntStep <- 1
    ret <- GRangesList()
    while (nrow(vcf_yield <- readVcf(tab, genome = genome, param = param))) {
        flog.info("Processing variants %i to %i...", cntVar + 1, cntVar + yieldSize)
        if (!(cntStep %% 10)) {
            flog.info("Position %s:%i", as.character(seqnames(vcf_yield)[1]), start(vcf_yield)[1])
        }
        mappingBias <- .calculateMappingBias(vcf_yield, min.normals, 
            min.normals.betafit, min.median.coverage.betafit, betafit.coverage.outliers)
        ret <- append(ret, GRangesList(mappingBias))
        cntVar <- cntVar + yieldSize
        cntStep <- cntStep + 1
    }
    bias <- .findMaxBetaShape(unlist(ret))
    attr(bias, "normal.panel.vcf.file") <- normal.panel.vcf.file
    attr(bias, "min.normals") <- min.normals
    attr(bias, "min.normals.betafit") <- min.normals.betafit
    attr(bias, "min.median.coverage.betafit") <- min.median.coverage.betafit
    attr(bias, "betafit.coverage.outliers") <- betafit.coverage.outliers
    attr(bias, "genome") <- genome
    bias
}
