# Make CMD check happy
globalVariables(names=c("..level.."))

#' Correct for GC bias
#' 
#' Takes as input coverage data and a mapping file for GC content. Will then
#' normalize coverage data for GC-bias.  Optionally plots the pre and
#' post normalization GC profiles.
#' 
#' 
#' @param coverage.file Coverage file or coverage data parsed with the
#' \code{\link{readCoverageFile}} function.
#' @param gc.gene.file File providing GC content for each exon in the coverage
#' files. First column in format CHR:START-END. Second column GC content (0 to
#' 1).  Third column provides gene symbols, which are optional, but used in
#' \code{\link{runAbsoluteCN}} to generate gene level calls. This file can be
#' generated with GATK GCContentByInterval tool or with the
#' \code{\link{calculateGCContentByInterval}} function.
#' @param output.file Optionally, write file with GC corrected coverage. Can be
#' read with the \code{\link{readCoverageFile}} function.
#' @param plot.gc.bias Optionally, plot GC profiles of the pre-normalized and
#' post-normalized coverage. Provides a quick visual check of coverage bias.
#' @param plot.max.density By default, if the number of intervals in the
#' probe-set is > 50000, uses a kernel density estimate to plot the coverage
#' distribution. This uses the \code{stat_density} function from the ggplot2
#' package. Using this parameter, change the threshold at which density
#' estimation is applied. If the \code{plot.gc.bias} parameter is set as
#' \code{FALSE}, this will be ignored.
#' @author Angad Singh, Markus Riester
#' @seealso \code{\link{calculateGCContentByInterval}}
#' @examples
#' 
#' normal.coverage.file <- system.file("extdata", "example_normal.txt", 
#'     package="PureCN")
#' gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
#'     package="PureCN")
#' coverage <- correctCoverageBias(normal.coverage.file, gc.gene.file)
#' 
#' @export correctCoverageBias
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line aes
#'             xlab ylab theme element_text facet_wrap stat_density2d
#'             scale_alpha_continuous scale_y_sqrt
#' @importFrom stats loess lm
#' @importFrom utils write.table
correctCoverageBias <- function(coverage.file, gc.gene.file,
output.file = NULL, plot.gc.bias = FALSE,
plot.max.density = 50000) {
    if (is.character(coverage.file)) {
        tumor  <- readCoverageFile(coverage.file)
    } else {
        tumor <- coverage.file
    }    
    
    tumor <- .addGCData(tumor, gc.gene.file, verbose=FALSE)

    ret <- .correctCoverageBiasLoess(tumor)

    coverage <- ret$coverage
    medDiploid <- ret$medDiploid
    #coverage <- .correctDuplicationBiasLoess(coverage)

    if (!is.null(output.file)) {
        .writeCoverage(coverage, output.file)
    }

    if (plot.gc.bias==TRUE) {
        if (length(coverage) < plot.max.density) {
            density <- "Low"
        } else {
            density <- "High"
        }
        tumor$norm_status <- "Pre-normalized"
        coverage$norm_status <- "Post-normalized"
        tumCov <- rbind(as.data.frame(tumor)[,c("coverage","average.coverage","gc_bias","norm_status")],
            as.data.frame(coverage)[,c("coverage","average.coverage","gc_bias","norm_status")])

        gcPlot <- tumCov[which(tumCov$average.coverage<quantile(tumCov$average.coverage,0.999, na.rm=TRUE)),]
        gcPlot$norm_status <- factor(gcPlot$norm_status, levels = c("Pre-normalized","Post-normalized"))
        plotMed <- medDiploid[,c("gcIndex","denom")]
        colnames(plotMed) <- c("gcIndex","gcNum")
        plotMed$norm_status <- "Pre-normalized"
        medDiploid$norm_status <- "Post-normalized"
        plotMed <- rbind(plotMed,medDiploid[,c("gcIndex","gcNum","norm_status")])
        plotMed$norm_status <- factor(plotMed$norm_status, 
            levels=c("Pre-normalized","Post-normalized"))

        if (density == "Low") {
            print(ggplot(gcPlot, aes_string(x="gc_bias", y="average.coverage")) + 
                geom_point(color='red', alpha=0.2) + 
                geom_line(data = plotMed, aes_string(x = 'gcIndex', y = 'gcNum'), color = 'blue') + 
                scale_y_sqrt() +
                xlab("GC content") + ylab("Coverage") + 
                theme(axis.text = element_text(size= 6), axis.title = element_text(size=16)) + 
                facet_wrap(~ norm_status, nrow=1))
        } else if (density == "High") {
            print(ggplot(gcPlot, aes_string(x="gc_bias", y="average.coverage")) + 
            geom_point(color="blue", alpha = 0.1) + 
            stat_density2d(aes(fill = ..level..), geom="polygon") + 
            scale_alpha_continuous(limits=c(0.1, 0), breaks=seq(0, 0.1, by = 0.025)) + 
            geom_line(data = plotMed, aes_string(x = 'gcIndex',y = 'gcNum'), color = 'red') + 
            scale_y_sqrt() +
            xlab("GC content") + ylab("Coverage") + 
            theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16)) + 
            facet_wrap(~norm_status, nrow=1))
        }
    }
    invisible(coverage)
}

#.correctDuplicationBiasLoess <- function(tumor) {
#    # ALPHA code
#    if (is.null(tumor$on.target)) tumor$on.target <- TRUE
#    if (is.null(tumor$duplication.rate)) return(tumor)
#    idx <- tumor$on.target
#    fit <- loess(tumor$duplication.rate[idx]~tumor$average.coverage[idx])
#    corFactor <- 1-predict(fit, tumor$average.coverage[idx])
#    tumor$coverage[idx] <- tumor$coverage[idx] / corFactor
#    tumor$counts[idx] <- tumor$counts[idx] / corFactor
#    .addAverageCoverage(tumor)
#}
    
.correctCoverageBiasLoess <- function(tumor) {
    if (is.null(tumor$on.target)) tumor$on.target <- TRUE
    gc_bias <- tumor$gc_bias
    for (on.target in c(FALSE, TRUE)) {
        tumor$valid <- tumor$on.target == on.target
        tumor$gc_bias <- gc_bias

        tumor$valid[tumor$average.coverage <= 0 | tumor$gc_bias < 0] <- FALSE

        if (!sum(tumor$valid)) next
        tumor$ideal <- TRUE
        routlier <- 0.01
        range <- quantile(tumor$average.coverage[tumor$valid], prob = 
            c(0, 1 - routlier), na.rm = TRUE)
        doutlier <- 0.001
        domain <- quantile(tumor$gc_bias[tumor$valid], prob = c(doutlier, 1 - doutlier), 
            na.rm = TRUE)
        
        tumor$ideal[!tumor$valid | 
            ( tumor$mappability < 1 & on.target ) |
            tumor$average.coverage <= range[1] |
            tumor$average.coverage > range[2] | 
            tumor$gc_bias < domain[1] | 
            tumor$gc_bias > domain[2]] <- FALSE
        
        if (!on.target) {    
            widthR <- quantile(width(tumor[tumor$ideal]), prob=0.1)
            tumor$ideal[width(tumor) < widthR] <- FALSE
        }
        rough <- loess(tumor$average.coverage[tumor$ideal] ~ tumor$gc_bias[tumor$ideal], 
            span = 0.03)
        i <- seq(0, 1, by = 0.001)
        final <- loess(predict(rough, i) ~ i, span = 0.3)
        cor.gc <- predict(final, tumor$gc_bias[tumor$valid])
        cor.gc.factor <- cor.gc/mean(tumor$average.coverage[tumor$ideal], na.rm=TRUE)
        cor.gc.factor[cor.gc.factor<=0] <- NA
        tumor$gc_bias <- as.integer(tumor$gc_bias*100)/100

        pre <- by(tumor$average.coverage[tumor$ideal], tumor$gc_bias[tumor$ideal], median, na.rm=TRUE)
        medDiploid <- as.data.frame(cbind(as.numeric(names(pre)),as.vector(pre)))
        colnames(medDiploid) <- c("gcIndex","denom")
        
#        tumor$average.coverage[tumor$valid] <- (tumor$average.coverage[tumor$valid] / cor.gc.factor)
        tumor$coverage[tumor$valid] <- (tumor$coverage[tumor$valid] / cor.gc.factor)
        tumor$counts[tumor$valid] <- (tumor$counts[tumor$valid] / cor.gc.factor)
        tumor <- .addAverageCoverage(tumor)

        post <- by(tumor$average.coverage[tumor$ideal], tumor$gc_bias[tumor$ideal], median, na.rm=TRUE)
        medDiploid$gcNum <- as.vector(post)
        tumor$ideal <- NULL
        tumor$valid <- NULL
    }
    ret <- list(coverage = tumor, medDiploid=medDiploid)
    ret
}

