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
#' @importFrom VGAM vglm Coef betabinomial dbetabinom
#' @export calculateMappingBiasVcf
calculateMappingBiasVcf <- function(normal.panel.vcf.file, min.normals = 2,
                                    min.normals.betafit = 7,
                                    min.median.coverage.betafit = 5,
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
            min.normals.betafit, min.median.coverage.betafit)
        ret <- append(ret, GRangesList(mappingBias))
        cntVar <- cntVar + yieldSize
        cntStep <- cntStep + 1
    }
    bias <- unlist(ret)
    attr(bias, "normal.panel.vcf.file") <- normal.panel.vcf.file
    attr(bias, "min.normals") <- min.normals
    attr(bias, "min.normals.betafit") <- min.normals.betafit
    attr(bias, "min.median.coverage.betafit") <- min.median.coverage.betafit
    attr(bias, "genome") <- genome
    bias
}

.calculateMappingBias <- function(nvcf, min.normals, min.normals.betafit = 7,
                                  min.median.coverage.betafit = 5) {
    if (ncol(nvcf) < 2) {
        .stopUserError("The normal.panel.vcf.file contains only a single sample.")
    }
    # TODO: deal with tri-allelic sites
    alt <- apply(geno(nvcf)$AD, c(1,2), function(x) x[[1]][2])
    ref <- apply(geno(nvcf)$AD, c(1,2), function(x) x[[1]][1])
    fa  <- apply(geno(nvcf)$AD, c(1,2), function(x) x[[1]][2]/sum(x[[1]]))
    x <- sapply(seq_len(nrow(nvcf)), function(i) {
        idx <- !is.na(fa[i,]) & fa[i,] > 0.05 & fa[i,] < 0.9
        shapes <- c(NA, NA)
        if (!sum(idx) >= min.normals) return(c(0, 0, 0, 0, shapes))
        dp <- alt[i,] + ref[i,] 
        mdp <- median(dp, na.rm = TRUE)   

        if (sum(idx) >= min.normals.betafit && mdp >= min.median.coverage.betafit) {
            fit <- suppressWarnings(try(vglm(cbind(alt[i,idx], 
                ref[i,idx]) ~ 1, betabinomial, trace = FALSE)))
            if (class(fit) == "try-error") {
                flog.warn("Could not fit beta binomial dist for %s (%s).",
                    as.character(rowRanges(nvcf[i])),
                    paste0(round(fa[i, idx], digits = 3), collapse=","))
            } else {    
                shapes <- Coef(fit)
            }
        }
        c(sum(ref[i,idx]), sum(alt[i,idx]), sum(idx), mean(fa[i, idx]), shapes)
    })
    # Add an average "normal" SNP (average coverage and allelic fraction > 0.4)
    # as empirical prior
    psMappingBias <- .adjustEmpBayes(x[1:4,]) * 2
    ponCntHits <- apply(geno(nvcf)$AD, 1, function(x)
        sum(!is.na(unlist(x))) / 2)

    tmp <- rowRanges(nvcf)
    mcols(tmp) <- NULL
    tmp$bias <- psMappingBias
    tmp$pon.count <- ponCntHits
    tmp$mu <- x[5,]
    tmp$rho <- x[6,]
    tmp
}
 
.readNormalPanelVcfLarge <- function(vcf, normal.panel.vcf.file,
    max.file.size=1, geno="AD", expand=FALSE) {
    genome <- genome(vcf)[1]
    if (!file.exists(normal.panel.vcf.file)) {
        .stopUserError("normal.panel.vcf.file ", normal.panel.vcf.file,
            " does not exist.")
    }
    if (file.size(normal.panel.vcf.file) / 1000^3 > max.file.size ||
        nrow(vcf) < 1000) {
        flog.info("Scanning %s...", normal.panel.vcf.file)
        nvcf <- readVcf(TabixFile(normal.panel.vcf.file), genome = genome,
            ScanVcfParam(which = rowRanges(vcf), info = NA, fixed = NA,
            geno = geno))
    } else {
        flog.info("Loading %s...", normal.panel.vcf.file)
        nvcf <- readVcf(normal.panel.vcf.file, genome = genome,
            ScanVcfParam(info = NA, fixed = NA, geno = geno))
        nvcf <- subsetByOverlaps(nvcf, rowRanges(vcf))
    }
    if (expand) nvcf <- expand(nvcf)
    nvcf
}

.adjustEmpBayes <- function(x) {
    # get all SNPs without dramatic bias
    xg <- x[, x[4, ] > 0.4, drop = FALSE]
    if (ncol(xg) < 2) {
        flog.warn("All SNPs in the database have significant mapping bias!%s",
            " Check your database.")
        shape1 <- 0
        shape2 <- 0
    } else {   
        # calculate the average number of ref and alt reads per sample
        shape1 <- sum(xg[1, ]) / sum(xg[3, ])
        shape2 <- sum(xg[2, ]) / sum(xg[3, ])
    }
    # add those as empirical bayes estimate to all SNPs
    x[1, ] <- x[1, ] + shape1
    x[2, ] <- x[2, ] + shape2
    # get the alt allelic fraction for all SNPs
    apply(x, 2, function(y) y[2] / sum(head(y,2)))
}
