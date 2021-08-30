#' Calculate Mapping Bias
#'
#' Function calculate mapping bias for each variant in the provided
#' panel of normals VCF.
#'
#'
#' @param normal.panel.vcf.file \code{character(1)} Combined VCF file of
#' a panel of normals, reference and alt counts as AD genotype field.
#' Needs to be compressed and indexed with bgzip and tabix, respectively.
#' @param min.normals Minimum number of normals with heterozygous SNP for
#' calculating position-specific mapping bias.
#' @param min.normals.betafit Minimum number of normals with heterozygous SNP
#' fitting a beta binomial distribution
#' @param min.normals.assign.betafit Minimum number of normals with
#' heterozygous SNPs to assign to a beta binomal fit cluster
#' @param min.median.coverage.betafit Minimum median coverage of normals with
#' heterozygous SNP for fitting a beta binomial distribution
#' @param num.betafit.clusters Maximum number of beta binomial fit clusters
#' @param min.betafit.rho Minimum dispersion factor rho
#' @param yieldSize See \code{TabixFile}
#' @param genome See \code{readVcf}
#' @return A \code{GRanges} object with mapping bias and number of normal
#' samples with this variant.
#' @author Markus Riester
#' @examples
#'
#' normal.panel.vcf <- system.file("extdata", "normalpanel.vcf.gz",
#'     package = "PureCN")
#' bias <- calculateMappingBiasVcf(normal.panel.vcf, genome = "h19")
#' saveRDS(bias, "mapping_bias.rds")
#'
#' @importFrom GenomicRanges GRangesList
#' @importFrom VGAM vglm Coef betabinomial dbetabinom
#' @importFrom data.table rbindlist
#' @importFrom mclust Mclust mclustBIC emControl
#' @export calculateMappingBiasVcf
calculateMappingBiasVcf <- function(normal.panel.vcf.file,
                                    min.normals = 1,
                                    min.normals.betafit = 7,
                                    min.normals.assign.betafit = 3,
                                    min.median.coverage.betafit = 5,
                                    num.betafit.clusters = 9,
                                    min.betafit.rho = 1e-04,
                                    yieldSize = 50000, genome) {
    tab <- TabixFile(normal.panel.vcf.file, yieldSize = yieldSize)
    open(tab)
    param <- ScanVcfParam(geno = c("AD"), fixed = "ALT", info = NA)
    cntVar <- 0
    cntStep <- 1
    ret <- GRangesList()
    while (nrow(vcf_yield <- readVcf(tab, genome = genome, param = param))) {
        flog.info("Processing variants %i to %i...", cntVar + 1, cntVar + yieldSize)
        if (!(cntStep %% 5)) {
            flog.info("Position %s:%i", as.character(seqnames(vcf_yield)[1]), start(vcf_yield)[1])
        }
        mappingBias <- .calculateMappingBias(nvcf = vcf_yield,
            min.normals = min.normals,
            min.normals.betafit = min.normals.betafit,
            min.normals.assign.betafit = min.normals.assign.betafit,
            min.median.coverage.betafit = min.median.coverage.betafit,
            num.betafit.clusters = num.betafit.clusters,
            min.betafit.rho = min.betafit.rho
        )
        if (length(ret)) {
            ret <- append(ret, GRangesList(mappingBias))
        } else {
            ret <- GRangesList(mappingBias)
        }
        cntVar <- cntVar + yieldSize
        cntStep <- cntStep + 1
    }
    bias <- unlist(ret)
    attr(bias, "normal.panel.vcf.file") <- normal.panel.vcf.file
    attr(bias, "min.normals") <- min.normals
    attr(bias, "min.normals.betafit") <- min.normals.betafit
    attr(bias, "min.normals.assign.betafit") <- min.normals.assign.betafit
    attr(bias, "min.median.coverage.betafit") <- min.median.coverage.betafit
    attr(bias, "num.betafit.clusters") <- num.betafit.clusters
    attr(bias, "min.betafit.rho") <- min.betafit.rho
    attr(bias, "genome") <- genome
    bias
}

#' Calculate Mapping Bias from GATK4 GenomicsDB
#'
#' Function calculate mapping bias for each variant in the provided
#' panel of normals GenomicsDB.
#'
#'
#' @param workspace Path to the GenomicsDB created by \code{GenomicsDBImport}
#' @param reference.genome Reference FASTA file.
#' @param min.normals Minimum number of normals with heterozygous SNP for
#' calculating position-specific mapping bias.
#' @param min.normals.betafit Minimum number of normals with heterozygous SNP
#' fitting a beta distribution
#' @param min.normals.assign.betafit Minimum number of normals with
#' heterozygous SNPs to assign to a beta binomal fit cluster
#' @param min.median.coverage.betafit Minimum median coverage of normals with
#' heterozygous SNP for fitting a beta distribution
#' @param num.betafit.clusters Maximum number of beta binomial fit clusters
#' @param min.betafit.rho Minimum dispersion factor rho
#' @param AF.info.field Field in the \code{workspace} that stores the allelic
#' fraction
#' @return A \code{GRanges} object with mapping bias and number of normal
#' samples with this variant.
#' @author Markus Riester
#' @examples
#'
#' \dontrun{
#' resources_file <- system.file("extdata", "gatk4_pon_db.tgz",
#'     package = "PureCN")
#' tmp_dir <- tempdir()
#' untar(resources_file, exdir = tmp_dir)
#' workspace <- file.path(tmp_dir, "gatk4_pon_db")
#' bias <- calculateMappingBiasGatk4(workspace, "hg19")
#' saveRDS(bias, "mapping_bias.rds")
#' unlink(tmp_dir, recursive=TRUE)
#' }
#'
#' @export calculateMappingBiasGatk4
#' @importFrom data.table dcast
#' @importFrom GenomeInfoDb rankSeqlevels
calculateMappingBiasGatk4 <- function(workspace, reference.genome,
                                    min.normals = 1,
                                    min.normals.betafit = 7,
                                    min.normals.assign.betafit = 3,
                                    min.median.coverage.betafit = 5,
                                    num.betafit.clusters = 9,
                                    min.betafit.rho = 1e-04,
                                    AF.info.field = "AF") {

    if (!requireNamespace("genomicsdb", quietly = TRUE) ||
        !requireNamespace("jsonlite", quietly = TRUE)
        ) {
        .stopUserError("Install the genomicsdb and jsonlite R packages for GenomicsDB import.")
    }
    workspace <- normalizePath(workspace, mustWork = TRUE)

    db <- genomicsdb::connect(workspace = workspace,
        vid_mapping_file = file.path(workspace, "vidmap.json"),
        callset_mapping_file = file.path(workspace, "callset.json"),
        reference_genome = reference.genome,
        c("DP", "AD", AF.info.field))

    jcallset <- jsonlite::read_json(file.path(workspace, "callset.json"))
    jvidmap <- jsonlite::read_json(file.path(workspace, "vidmap.json"))
    
    # get all available arrays
    arrays <- sapply(dir(workspace, full.names = TRUE), file.path, "genomicsdb_meta_dir")
    arrays <- basename(names(arrays)[which(file.exists(arrays))])
    # get offsets and lengths
    contigs <- sapply(arrays, function(ary) strsplit(ary, "\\$")[[1]][1])
    contigs <- jvidmap$contigs[match(contigs, sapply(jvidmap$contigs, function(x) x$name))]
    idx <- order(rankSeqlevels(sapply(contigs, function(x) x$name)))
    row_ranges <- list(range(sapply(jcallset$callsets, function(x) x$row_idx)))
    flog.info("Found %i contigs and %i columns.",
        length(contigs), row_ranges[[1]][2] - row_ranges[[1]][1] + 1)

    parsed_ad_list <- lapply(idx, function(i) {
        c_offset <- as.numeric(contigs[[i]]$tiledb_column_offset)
        c_length <- as.numeric(contigs[[i]]$length)

        flog.info("Processing %s (offset %.0f, length %.0f)...",
            arrays[i], c_offset, c_length)
        query <- data.table(genomicsdb::query_variant_calls(db,
            array = arrays[i],
            column_ranges = list(c(c_offset, c_offset + c_length)),
            row_ranges = row_ranges
        ))
        .parseADGenomicsDb(query, AF.info.field)
    })
    genomicsdb::disconnect(db)
    flog.info("Collecting variant information...")
    # concat with minimal memory overhead
    m_alt <- as.matrix(rbindlist(lapply(parsed_ad_list, function(x) data.frame(x$alt)), fill = TRUE))
    for (i in seq_along(parsed_ad_list)) parsed_ad_list[[i]]$alt <- NULL
    m_ref <- as.matrix(rbindlist(lapply(parsed_ad_list, function(x) data.frame(x$ref)), fill = TRUE))
    for (i in seq_along(parsed_ad_list)) parsed_ad_list[[i]]$ref <- NULL
    gr <- unlist(GRangesList(lapply(parsed_ad_list, function(x) x$gr)))
    parsed_ad_list <- NULL
    bias <- .calculateMappingBias(nvcf = NULL,
        alt = m_alt, ref = m_ref, gr = gr,
        min.normals = min.normals,
        min.normals.betafit = min.normals.betafit,
        min.normals.assign.betafit = min.normals.assign.betafit,
        min.median.coverage.betafit = min.median.coverage.betafit,
        num.betafit.clusters = num.betafit.clusters,
        min.betafit.rho = min.betafit.rho,
        verbose = TRUE
    )
    attr(bias, "workspace") <- workspace
    attr(bias, "min.normals") <- min.normals
    attr(bias, "min.normals.betafit") <- min.normals.betafit
    attr(bias, "min.median.coverage.betafit") <- min.median.coverage.betafit
    return(bias)
}

.parseADGenomicsDb <- function(query, AF.info.field = "AF") {
    ref <-  dcast(query, CHROM + POS + END + REF + ALT ~ SAMPLE, value.var = "AD")
    af <-  dcast(query, CHROM + POS + END + REF + ALT ~ SAMPLE, value.var = AF.info.field)
    gr <- GRanges(seqnames = ref$CHROM, IRanges(start = ref$POS, end = ref$END),
        strand = NULL, DataFrame(REF = ref$REF, ALT = ref$ALT))
    genomic_change <- paste0(as.character(gr), "_", ref$REF, ">", ref$ALT)
    ref <- as.matrix(ref[,-(1:5)])
    af <- as.matrix(af[,-(1:5)])
    alt <- round(ref / (1 - af) - ref)
    rownames(ref) <- genomic_change
    rownames(af) <- genomic_change
    rownames(alt) <- genomic_change
    list(ref = ref, alt = alt, gr = gr)
}
    
.calculateMappingBias <- function(nvcf, alt = NULL, ref = NULL, gr = NULL,
                                  min.normals, min.normals.betafit = 7,
                                  min.normals.assign.betafit = 3,
                                  min.median.coverage.betafit = 5,
                                  num.betafit.clusters = 9,
                                  min.betafit.rho = 1e-04,
                                  verbose = FALSE) {
    if (!is.null(nvcf)) {
        if (ncol(nvcf) < 2) {
            .stopUserError("The normal.panel.vcf.file contains only a single sample.")
        }
        # TODO: deal with tri-allelic sites
        alt <- apply(geno(nvcf)$AD, c(1, 2), function(x) x[[1]][2])
        ref <- apply(geno(nvcf)$AD, c(1, 2), function(x) x[[1]][1])
        fa  <- apply(geno(nvcf)$AD, c(1, 2), function(x) x[[1]][2] / sum(x[[1]]))
        gr <- rowRanges(nvcf)
    } else {
        if (is.null(alt) || is.null(ref) || is.null(gr) || ncol(ref) != ncol(alt)) {
            .stopRuntimeError("Either nvcf or valid alt and ref required.")
        }
        if (ncol(alt) < 2) {
            .stopUserError("The normal.panel.vcf.file contains only a single sample.")
        }
        fa <- alt / (ref + alt)
    }
    ponCntHits <- apply(alt, 1, function(x) sum(!is.na(x)))
    if (verbose) {
        flog.info("Fitting beta-binomial distributions. Might take a while...")
    }
    x <- sapply(seq_len(nrow(fa)), function(i) {
        if (verbose && !(i %% 50000)) {
            flog.info("Position %s (variant %.0f/%.0f)", as.character(gr)[i], i, nrow(alt))
        }
        idx <- !is.na(fa[i, ]) & fa[i, ] > 0.05 & fa[i, ] < 0.9
        shapes <- c(NA, NA)
        if (!sum(idx) >= min.normals) return(c(0, 0, 0, 0, shapes))
        dp <- alt[i, ] + ref[i, ]
        mdp <- median(dp, na.rm = TRUE)

        if (sum(idx) >= min.normals.betafit && mdp >= min.median.coverage.betafit) {
            fit <- suppressWarnings(try(vglm(cbind(alt[i, idx],
                ref[i, idx]) ~ 1, betabinomial, trace = FALSE)))
            if (is(fit, "try-error")) {
                flog.warn("Could not fit beta binomial dist for %s (%s).",
                    as.character(gr)[i],
                    paste0(round(fa[i, idx], digits = 3), collapse = ","))
            } else {
                shapes <- Coef(fit)
            }
        }
        c(sum(ref[i, idx]), sum(alt[i, idx]), sum(idx), mean(fa[i, idx]), shapes)
    })
    # Add an average "normal" SNP (average coverage and allelic fraction > 0.4)
    # as empirical prior
    gr$bias <- .adjustEmpBayes(x[1:4, ]) * 2
    gr$pon.count <- ponCntHits
    gr$mu <- x[5, ]
    gr$rho <- pmax(min.betafit.rho, x[6, ])
    gr <- .clusterFa(gr, fa, alt, ref, min.normals.assign.betafit, num.betafit.clusters)
    gr <- gr[order(gr$pon.count, decreasing = TRUE)]
    gr <- sort(gr)
    gr$triallelic <- FALSE
    gr$triallelic[duplicated(gr, fromLast = FALSE) |
                  duplicated(gr, fromLast = TRUE)] <- TRUE
    gr
}
 
.readNormalPanelVcfLarge <- function(vcf, normal.panel.vcf.file,
    max.file.size = 1, geno = "AD", expand = FALSE) {
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
    apply(x, 2, function(y) y[2] / sum(head(y, 2)))
}

.clusterFa <- function(gr, fa, alt, ref, min.normals.assign.betafit = 3, num.betafit.clusters = 9) {
    idx <- !is.na(gr$mu)
    gr$clustered <- FALSE
    if (sum(idx) < num.betafit.clusters) return(gr)
    flog.info("Clustering beta binomial fits...")
    fit <- Mclust(mcols(gr)[idx, c("mu", "rho")], G = seq_len(num.betafit.clusters))
    x <- sapply(seq_len(nrow(fa)), function(i) {
        idx <- !is.na(fa[i, ]) & fa[i, ] > 0.05 & fa[i, ] < 0.9
        if (!sum(idx) >= min.normals.assign.betafit) return(NA)
        ll <- apply(fit$parameters$mean, 2, function(m) {
                sum(dbetabinom(x = alt[i, idx], prob = m[1],
                   size = alt[i, idx] + ref[i, idx],
                   rho = m[2], log = TRUE))
        })
        if (is.infinite(max(ll))) return(NA)
        which.max(ll)
    })
    idxa <- is.na(gr$mu) & !is.na(x)
    flog.info("Assigning (%i/%i) variants a clustered beta binomal fit.",
        sum(idxa), length(gr))
    gr$mu[idxa] <- fit$parameters$mean[1, x[idxa]]
    gr$rho[idxa] <- fit$parameters$mean[2, x[idxa]]
    gr$clustered[idxa] <- TRUE
    return(gr)
}
