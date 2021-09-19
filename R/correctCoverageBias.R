# Make CMD check happy
globalVariables(names = c("..level.."))

#' Correct for library-specific coverage biases
#'
#' Takes as input coverage data and a mapping file for GC content and
#' optionally replication timing. Will then normalize coverage data for
#' GC-bias.  Plots the pre and post normalization GC profiles.
#'
#'
#' @param coverage.file Coverage file or coverage data parsed with the
#' \code{\link{readCoverageFile}} function.
#' @param interval.file File providing GC content for each exon in the coverage
#' files. First column in format CHR:START-END. Additional optional columns
#' provide gene symbols, mappability and replication timing. This file is
#' generated with the \code{\link{preprocessIntervals}} function.
#' @param output.file Optionally, write file with GC corrected coverage. Can be
#' read with the \code{\link{readCoverageFile}} function.
#' @param plot.bias Optionally, plot profiles of the pre-normalized and
#' post-normalized coverage. Provides a quick visual check of coverage bias.
#' @param plot.max.density By default, if the number of intervals in the
#' probe-set is > 50000, uses a kernel density estimate to plot the coverage
#' distribution. This uses the \code{stat_density} function from the ggplot2
#' package. Using this parameter, change the threshold at which density
#' estimation is applied. If the \code{plot.bias} parameter is set as
#' \code{FALSE}, this will be ignored.
#' @param output.qc.file Write miscellaneous coverage QC metrics to file.
#' @author Angad Singh, Markus Riester
#' @seealso \code{\link{preprocessIntervals}}
#' @examples
#'
#' normal.coverage.file <- system.file("extdata", "example_normal.txt",
#'     package = "PureCN")
#' interval.file <- system.file("extdata", "example_intervals.txt",
#'     package = "PureCN")
#' coverage <- correctCoverageBias(normal.coverage.file, interval.file)
#'
#' @export correctCoverageBias
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line aes
#'             xlab ylab theme element_text facet_wrap stat_density2d
#'             scale_alpha_continuous scale_y_sqrt geom_abline
#'             coord_trans
#' @importFrom gridExtra grid.arrange
#' @importFrom stats loess lm predict
#' @importFrom data.table fwrite
correctCoverageBias <- function(coverage.file, interval.file,
output.file = NULL, plot.bias = FALSE, plot.max.density = 50000,
output.qc.file = NULL) {

    if (is.character(coverage.file)) {
        raw  <- readCoverageFile(coverage.file)
    } else {
        raw <- coverage.file
    }

    if (max(raw$average.coverage[raw$on.target], na.rm = TRUE) <= 0) {
        .stopUserError("Provided coverage is zero, most likely due to a corrupt BAM file.")
    }

    raw <- .addGCData(raw, interval.file, verbose = FALSE)
    ret <- .correctCoverageBiasLoess(raw)
    if (plot.bias) {
        gp1 <- .plotGcBias(raw, ret$coverage, ret$medDiploid, plot.max.density)
    }

    gc <- ret$coverage
    ret <- .correctRepTimingBiasLinear(gc)
    if (plot.bias) {
        # reptiming available?
        if (!is.null(ret$lmFit)) {
            gp2 <- .plotRepBias(gc, ret$coverage, ret$lmFit, plot.max.density)
            grid.arrange(gp1, gp2, nrow = 2)
        } else {
            print(gp1)
        }
    }
    if (!is.null(output.file)) {
        .writeCoverage(ret$coverage, output.file)
    }
    if (!is.null(output.qc.file)) {
        .writeQCFile(raw, gc, ret$coverage, output.qc.file)
    }
    invisible(ret$coverage)
}

.MoM <- function(x, plot = FALSE) {
    if (length(x) < 10) return(NA)
    x <- median(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
}

.writeQCFile <- function(raw, gc, final, output.qc.file) {
    mom <- unlist(lapply(list(raw, gc, final), function(x) sapply(c(TRUE, FALSE), function(b)
        .MoM(x[which(x$on.target == b)]$average.coverage))))

    meanOn     <- mean(raw[which(raw$on.target)]$average.coverage, na.rm = TRUE)
    meanOff    <- mean(raw[which(!raw$on.target)]$average.coverage, na.rm = TRUE)
    meanDupOn  <- mean(raw[which(raw$on.target)]$duplication.rate, na.rm = TRUE)
    meanDupOff <- mean(raw[which(!raw$on.target)]$duplication.rate, na.rm = TRUE)

    qc <- c(meanOn, meanOff, meanDupOn, meanDupOff, mom)

    qc <- data.frame(matrix(qc, nrow = 1))
    colnames(qc)[1:2] <- c("mean.coverage.ontarget", "mean.coverage.offtarget")
    colnames(qc)[3:4] <- c("mean.duplication.ontarget", "mean.duplication.offtarget")
    colnames(qc)[5:10] <- paste0("mom.", c("raw", "raw", "post.gc", "post.gc",
                             "post.reptiming", "post.reptiming"),
                             ".", rep(c("ontarget", "offtarget"), 3))
    fwrite(qc, file = output.qc.file, row.names = FALSE, quote = FALSE,
        sep = " ")
}
    
.createCoverageGgplot <- function(raw, normalized, plot.max.density, x, log = FALSE) {
    if (length(normalized) < plot.max.density) {
        density <- "Low"
    } else {
        density <- "High"
    }
    raw$norm_status <- "Pre-normalized"
    normalized$norm_status <- "Post-normalized"
    ids <- c("coverage", "average.coverage", "gc_bias", "reptiming",
             "norm_status", "on.target")
    tumCov <- rbind(as.data.frame(raw)[, ids],
                    as.data.frame(normalized)[, ids])

    tumCov <- tumCov[which(tumCov$average.coverage <
        quantile(tumCov$average.coverage, 0.999, na.rm = TRUE)), ]

    if (sum(!tumCov$on.target)) {
        tumCov <- tumCov[which(tumCov$on.target | tumCov$average.coverage <
            quantile(tumCov$average.coverage[!tumCov$on.target], 0.99, na.rm = TRUE)), ]
    }
    tumCov$on.target <- factor(ifelse(tumCov$on.target, "on-target", "off-target"),
        levels = c("on-target", "off-target"))
    tumCov$norm_status <- factor(tumCov$norm_status,
        levels = c("Pre-normalized", "Post-normalized"))
    if (log) tumCov$average.coverage <- log(tumCov$average.coverage)

    gp <- ggplot(tumCov, aes_string(x = x, y = "average.coverage"))
    if (density == "Low") {
        gp <- gp + geom_point(color = "red", alpha = 0.2)
    } else if (density == "High") {
        gp <- gp + geom_point(color = "blue", alpha = 0.1) +
            stat_density2d(aes(fill = ..level..), geom = "polygon") +
            scale_alpha_continuous(limits = c(0.1, 0),
                                   breaks = seq(0, 0.1, by = 0.025))
    }

    gp + ylab(paste0(if (log) "Log-" else "", "Coverage")) +
       facet_wrap(~on.target + norm_status, ncol = 2, scales = "free_y")
}

.plotGcBias <- function(raw, normalized, medDiploid, plot.max.density) {
    gp <- .createCoverageGgplot(raw, normalized, plot.max.density, "gc_bias")

    plotMed <- medDiploid[, c("gcIndex", "denom")]
    colnames(plotMed) <- c("gcIndex", "gcNum")
    plotMed$norm_status <- "Pre-normalized"
    medDiploid$norm_status <- "Post-normalized"
    plotMed <- rbind(plotMed, medDiploid[, c("gcIndex", "gcNum", "norm_status")])
    plotMed$norm_status <- factor(plotMed$norm_status,
        levels = c("Pre-normalized", "Post-normalized"))
    plotMed$on.target <- factor("on-target", levels = c("on-target", "off-target"))
    gp <- gp + geom_line(data = plotMed, aes_string(x = "gcIndex", y = "gcNum"), color = "blue") +
          xlab("GC content")
    gp
}

.plotRepBias <- function(raw, normalized, lmFit, plot.max.density) {
    gp <- .createCoverageGgplot(raw, normalized, plot.max.density, "reptiming", log = TRUE) +
            xlab("Replication Timing")
    plotMed <- do.call(rbind, lapply(lmFit, function(x)
        data.frame(rbind(x$before$coefficients, x$after$coefficients),
            norm_status = c("Pre-normalized", "Post-normalized"),
            on.target = x$on.target)))
    colnames(plotMed)[1:2] <- c("intercept", "slope")
    plotMed$norm_status <- factor(plotMed$norm_status,
        levels = c("Pre-normalized", "Post-normalized"))
    plotMed$on.target <- factor(ifelse(plotMed$on.target, "on-target", "off-target"),
        levels = c("on-target", "off-target"))
    gp <- gp + geom_abline(data = plotMed,
        aes_string(intercept = "intercept", slope = "slope"), color = "blue")
    gp
}

.correctRepTimingBiasLinear <- function(tumor) {
    # ALPHA code
    if (is.null(tumor$on.target)) tumor$on.target <- TRUE
    if (is.null(tumor$reptiming) || sum(!is.na(tumor$reptiming)) < 100) {
        return(list(coverage = tumor, lmFit = NULL))
    }
    doutlier <- 0.001
    domain <- quantile(tumor$reptiming, probs = c(doutlier, 1 - doutlier), na.rm = TRUE)

    lmFit <- list()

    for (on.target in c(FALSE, TRUE)) {
        # ignore problematic intervals (missing or outlier data)
        idx <- tumor$on.target == on.target &
               !is.na(tumor$reptiming) &
               !is.na(tumor$average.coverage) &
               tumor$average.coverage > 0 &
               tumor$reptiming >= domain[1] &
               tumor$reptiming <= domain[2]
        if (!sum(idx)) next
        fit <- lm(log(tumor$average.coverage[idx])~tumor$reptiming[idx])
        corFactor <- exp(predict(fit) - mean(predict(fit)))

        tumor$coverage[idx] <- tumor$coverage[idx] / corFactor
        tumor$counts[idx] <- tumor$counts[idx] / corFactor
        tumor <- .addAverageCoverage(tumor)
        fitAfter <- lm(log(tumor$average.coverage[idx])~tumor$reptiming[idx])
        lmFit[[length(lmFit) + 1]] <- list(before = fit,
            after = fitAfter, on.target = on.target)
    }
    list(coverage = tumor, lmFit = lmFit)
}
    
.correctCoverageBiasLoess <- function(tumor) {
    if (is.null(tumor$on.target)) tumor$on.target <- TRUE
    gc_bias <- tumor$gc_bias
    for (on.target in c(FALSE, TRUE)) {
        tumor$valid <- tumor$on.target == on.target
        tumor$gc_bias <- gc_bias

        tumor$valid[tumor$average.coverage <= 0 | tumor$gc_bias < 0] <- FALSE

        if (!any(tumor$valid)) next
        tumor$ideal <- TRUE
        routlier <- 0.01
        range <- quantile(tumor$average.coverage[tumor$valid], prob =
            c(0, 1 - routlier), na.rm = TRUE)
        doutlier <- 0.001
        domain <- quantile(tumor$gc_bias[tumor$valid], prob = c(doutlier, 1 - doutlier),
            na.rm = TRUE)

        tumor$ideal[!tumor$valid |
            (tumor$mappability < 1 & on.target) |
            tumor$average.coverage <= range[1] |
            tumor$average.coverage > range[2] |
            tumor$gc_bias < domain[1] |
            tumor$gc_bias > domain[2]] <- FALSE

        if (!any(tumor$ideal)) next

        if (!on.target) {
            widthR <- quantile(width(tumor[tumor$ideal]), prob = 0.1)
            tumor$ideal[width(tumor) < widthR] <- FALSE
        }
        rough <- loess(tumor$average.coverage[tumor$ideal] ~ tumor$gc_bias[tumor$ideal],
            span = 0.03)
        i <- seq(0, 1, by = 0.001)
        final <- loess(predict(rough, i) ~ i, span = 0.3)
        cor.gc <- predict(final, tumor$gc_bias[tumor$valid])
        cor.gc.factor <- cor.gc / mean(tumor$average.coverage[tumor$ideal], na.rm = TRUE)
        cor.gc.factor[cor.gc.factor <= 0] <- NA
        tumor$gc_bias <- as.integer(tumor$gc_bias * 100) / 100

        pre <- by(tumor$average.coverage[tumor$ideal], tumor$gc_bias[tumor$ideal], median, na.rm = TRUE)
        medDiploid <- as.data.frame(cbind(as.numeric(names(pre)), as.vector(pre)))
        colnames(medDiploid) <- c("gcIndex", "denom")

        tumor$coverage[tumor$valid] <- (tumor$coverage[tumor$valid] / cor.gc.factor)
        tumor$counts[tumor$valid] <- (tumor$counts[tumor$valid] / cor.gc.factor)
        tumor <- .addAverageCoverage(tumor)

        post <- by(tumor$average.coverage[tumor$ideal], tumor$gc_bias[tumor$ideal], median, na.rm = TRUE)
        medDiploid$gcNum <- as.vector(post)
        tumor$ideal <- NULL
        tumor$valid <- NULL
    }
    list(coverage = tumor, medDiploid = medDiploid)
}
