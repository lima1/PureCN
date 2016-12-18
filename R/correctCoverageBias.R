# Make CMD check happy
globalVariables(names=c("gcIndex", "gcNum", "..level.."))

correctCoverageBias <- structure(function(# Correct for GC bias
### Takes as input coverage data in GATK format (or data 
### read by \code{\link{readCoverageGatk}}) and a mapping file for 
### GC content, and normalize coverage data for bias correction. 
### Optionally plots the pre and post normalization GC profiles.
coverage.file, 
### Exon coverage file as produced by GATK. Either a file name
### or data parsed with the \code{\link{readCoverageGatk}} function.
gc.gene.file,
### File providing GC content for each exon in the coverage files.
### First column in format CHR:START-END. Second column GC content (0 to 1). 
### Third column provides gene symbols, which are optional, but used in
### \code{\link{runAbsoluteCN}} to generate gene level calls. This file
### can be generated with GATK GCContentByInterval tool or with the
### \code{\link{calculateGCContentByInterval}} function.
##seealso<< \code{\link{calculateGCContentByInterval}}
output.file=NULL,
### Optionally, write file with GC corrected coverage. Can be read with
### the \code{\link{readCoverageGatk}} function.
method=c("LOESS","POLYNOMIAL"),
### Two options for normalization are available: The default "LOESS" largely
### follows the GC correction of the TitanCNA package. The "POLYNOMIAL" method
### models the coverage data as a polynomial of degree three and normalizes
### using the EM approach. The "POLYNOMIAL" is expected to be more robust for
### smaller targeted panels.
plot.gc.bias=FALSE,
### Optionally, plot GC profiles of the pre-normalized and post-normalized
### coverage. Provides a quick visual check of coverage bias.
plot.max.density=50000
### By default, if the number of intervals in the probe-set is > 50000, uses
### a kernel density estimate to plot the coverage distribution. This uses
### the \code{stat_density} function from the ggplot2 package. Using this
### parameter, change the threshold at which density estimation is applied.
### If the \code{plot.gc.bias} parameter is set as \code{FALSE}, this will be
### ignored.
) {
    if (is.character(coverage.file)) {
        tumor  <- readCoverageGatk(coverage.file)
    } else {
        tumor <- coverage.file
    }    
    
    gc <- .loadGcGeneFile(gc.gene.file, tumor)

    method <- match.arg(method)
    if (method=="LOESS") {
        ret <- .correctCoverageBiasLoess(tumor, gc)
    } else if(method=="POLYNOMIAL") {
        ret <- .correctCoverageBiasPolynomial(tumor, gc)
    } else {
        .stopUserError(
        "Invalid method specified. Please specify one of ",
        "\"LOESS\" or \"POLYNOMIAL\"."
        )
    }
    coverage <- ret$coverage
    medDiploid <- ret$medDiploid

    if (!is.null(output.file)) {
        tmp <- coverage[, c("probe", "coverage", "average.coverage")]
        colnames(tmp) <- c("Target", "total_coverage", "average_coverage")
        write.table(tmp, file=output.file, row.names=FALSE, quote=FALSE)
    }

    if (plot.gc.bias==TRUE) {
        if (nrow(gc) < plot.max.density) {
            density <- "Low"
        } else {
            density <- "High"
        }
        tumor$gc_bias <- gc$gc_bias
        tumor$norm_status <- "Pre-normalized"
        coverage$gc_bias <- gc$gc_bias
        coverage$norm_status <- "Post-normalized"
        tumCov <- rbind(tumor[,c("probe","coverage","average.coverage","gc_bias","norm_status")],
            coverage[,c("probe","coverage","average.coverage","gc_bias","norm_status")])

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
                geom_line(data = plotMed, aes(x = gcIndex, y = gcNum), color = 'blue') + 
                xlab("GC content") + ylab("Coverage") + 
                theme(axis.text = element_text(size= 6), axis.title = element_text(size=16)) + 
                facet_wrap(~ norm_status, nrow=1))
        } else if (density == "High") {
            print(ggplot(gcPlot, aes_string(x="gc_bias", y="average.coverage")) + 
            geom_point(color="blue", alpha = 0.1) + 
            stat_density2d(aes(fill = ..level..), geom="polygon") + 
            scale_alpha_continuous(limits=c(0.1, 0), breaks=seq(0, 0.1, by = 0.025)) + 
            geom_line(data = plotMed, aes(x = gcIndex,y = gcNum), color = 'red') + 
            xlab("GC content") + ylab("Coverage") + 
            theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16)) + 
            facet_wrap(~norm_status, nrow=1))
        }
    }
##author<< Angad Singh,
    invisible(coverage)
### GC normalized coverage.
}, ex=function() {
normal.coverage.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
    package="PureCN")
# normalize using default LOESS method
coverage <- correctCoverageBias(normal.coverage.file, gc.gene.file)
# normalize with POLYNOMIAL method for small panels
coverage <- correctCoverageBias(normal.coverage.file, gc.gene.file, 
    method="POLYNOMIAL", plot.gc.bias=TRUE)
})

.loadGcGeneFile <- function(gc.gene.file, tumor) {
    gc <- read.delim(gc.gene.file)

    if (is.null(gc$gc_bias)) {
        .stopUserError("gc.gene.file header invalid.")
    }
    
    if (!identical(as.character(gc[,1]), as.character(tumor[,1]))) {
        if (sum(!as.character(tumor[,1]) %in% as.character(gc[,1])) > 0) {
            .stopUserError(
            "Intervals in coverage.file and gc.gene.file different.\n",
            "Some intervals in coverage have no GC information.\n",
            "Please re-calculate coverage using the intervals specified ",
            "in gc.gene.file or provide the correct gc.gene.file."
            )
        }
        warning(
        "Interval files in coverage.file and gc.gene.file different.")
        gc <- gc[match(as.character(tumor[,1]), as.character(gc[,1])),]
    }
    gc
}

.correctCoverageBiasLoess <- function(tumor, gc) {
    # taken from TitanCNA
    gc$valid <- TRUE
    gc$valid[tumor$average.coverage <= 0 | gc$gc_bias < 0] <- FALSE
    gc$ideal <- TRUE
    routlier <- 0.01
    range <- quantile(tumor$average.coverage[gc$valid], prob = 
        c(0, 1 - routlier), na.rm = TRUE)
    doutlier <- 0.001
    domain <- quantile(gc$gc_bias[gc$valid], prob = c(doutlier, 1 - doutlier), 
        na.rm = TRUE)
    gc$ideal[!gc$valid | 
        tumor$average.coverage <= range[1] |
        tumor$average.coverage > range[2] | 
        gc$gc_bias < domain[1] | 
        gc$gc_bias > domain[2]] <- FALSE

    rough <- loess(tumor$average.coverage[gc$ideal] ~ gc$gc_bias[gc$ideal], 
        span = 0.03)
    i <- seq(0, 1, by = 0.001)
    final <- loess(predict(rough, i) ~ i, span = 0.3)
    cor.gc <- predict(final, gc$gc_bias)
    cor.gc.factor <- cor.gc/mean(tumor$average.coverage, na.rm=TRUE)

    tumor$gc_bias <- as.integer(gc$gc_bias*100)/100
    pre <- by(tumor$average.coverage[gc$ideal], tumor$gc_bias[gc$ideal], median, na.rm=TRUE)
    medDiploid <- as.data.frame(cbind(as.numeric(names(pre)),as.vector(pre)))
    colnames(medDiploid) <- c("gcIndex","denom")

    tumor$average.coverage <- tumor$average.coverage / cor.gc.factor
    tumor$coverage <- tumor$coverage / cor.gc.factor

    post <- by(tumor$average.coverage[gc$ideal], tumor$gc_bias[gc$ideal], median, na.rm=TRUE)
    medDiploid$gcNum <- as.vector(post)
    tumor$gc_bias <- NULL

    ret <- list(coverage = tumor, medDiploid=medDiploid)
    ret
}

.correctCoverageBiasPolynomial <- function(tumor, gc) {
    coverage <- .cleanupSample(tumor,gc)
    coveragePloidy <- .inferPloidyRegions(coverage)
    coverageNormalized <- .normalizeRunningMedian(coveragePloidy)
    coverageNormalized
}

.cleanupSample <- function(tumor, gc) {
  if(dim(tumor)[1]==0|dim(gc)[1]==0|!identical(as.character(gc[,1]), as.character(tumor[,1]))) {
    .stopUserError(
    paste("Mismatching interval sets:",dim(tumor),dim(gc),sep=" ")
    )
  }
min_cov <- 10
  tumor$length <- tumor$probe_end - tumor$probe_start + 1
  coverage <- merge(x=gc,y=tumor,by.x="Target",by.y="probe",sort=FALSE)[,
    c("Target","gc_bias","average.coverage","length")]

  # Remove instances of coverage < min_cov in sample or zero coverage, as well as chrX/Y
  coverage$usable <- TRUE
  coverage[is.na(coverage$average.coverage),"usable"] <- FALSE
  coverage[grep('chrY',as.character(coverage$Target)),"usable"] <- FALSE
  coverage[which(coverage$average.coverage<min_cov),"usable"] <- FALSE
  coverage[which(coverage$regLength<75),"usable"] <- FALSE

  coverage
}

.inferPloidyRegions <- function (coverage) {
  x <- coverage$gc[which(coverage$usable==TRUE)]
  y <- coverage$average.coverage[which(coverage$usable==TRUE)]

  # Need to init polynomial parameters here. Instead of randomly initializing, run
  # regression on ALL data points, not excluding on basis of assigned copy numbers or others.
  reg <- lm(y ~ poly(x, 3, raw=TRUE))
  ploidy <- 2
  P <- ploidy
  oldp0 <- 0 # "old" polynomial parameters retained to check for convergence
  oldp1 <- 0
  oldp2 <- 0
  oldp3 <- 0
  p0 <- summary(reg)$coefficients[1,1]
  p1 <- summary(reg)$coefficients[2,1]
  p2 <- summary(reg)$coefficients[3,1]
  p3 <- summary(reg)$coefficients[4,1]
  oldsum <- 0
  itercount=1
  # Run loop until parameters converge
  while(1) {
    newsum <- sum(p3*x*x*x+p2*x*x+p1*x+p0)+sum(oldp3*x*x*x+oldp2*x*x+oldp1*x+oldp0)
    if(abs(newsum-oldsum)<5) break
        itercount=itercount+1
        if(itercount==100) break

    oldsum <- newsum
    oldp0 <- p0
    oldp1 <- p1
    oldp2 <- p2
    oldp3 <- p3

    # Assign copy numbers to all data points: k' = argmin(k) abs(f(g_i)*k/P - RC_i)
max_cnv <- 10
    cn <- array(data=NA,c(max_cnv+1,length(x)))
    for(k in 0:max_cnv) {
      f_xi <- p3*x*x*x+p2*x*x+p1*x+p0
      cn[k+1,] <- abs(f_xi*k/P - y)
    }
    # Assign CNVs in current iteration
    min_index <- lapply(1:ncol(cn), function(i) which.min(cn[,i])-1)
    coverage[which(coverage$usable==TRUE),"CNV"] <- sapply(min_index, function(x){as.numeric(x[1])})

    medx <- seq(0.25, 0.80, 0.05) # Full GC range for fitting polynomial
    medy <- rep(NA, length(medx)) # Median coverage in given GC window
    if(length(x)<10000) medi <- seq(2,9) #Work on GC range (0.30,0.65) only. Beyond that all points are included.
    if(length(x)>=10000) medi <- seq(2,11) #Work on GC range (0.30,0.75) only. Beyond that all points are included.

    # Select only points with given ploidy in (0.30,0.65) and all points outside that. selects 2 classes of points:
        # regions which have the expected ploidy (diploid or user provided sample ploidy)
        # regions with G+C content below 0.3 or above 0.65
    topp <- quantile(coverage$average.coverage,0.99,na.rm=TRUE)
    botp <- quantile(coverage$average.coverage,0.01,na.rm=TRUE)
    coverage[which(coverage$usable==TRUE),"diploid"] <- FALSE
    coverage[which(coverage$CNV==P & coverage$gc_bias>=min(medx[medi]) & coverage$gc_bias<max(medx[medi]) & coverage$average.coverage<topp & coverage$average.coverage>botp & coverage$usable==TRUE),"diploid"] <- TRUE

    medcount <- rep(NA,length(medx)) # #points in window, used to weight the regression equation
    for(i in seq(1,length(medx))) {
      gc <- medx[i] # GC window of 0.05
      gc_min <- gc-.025
      gc_max <- gc+.025
      medgc <- coverage[which(coverage$gc_bias>=gc_min & coverage$gc_bias<gc_max & coverage$diploid==TRUE),] # All points in given GC window
      medy[i] <- median(medgc$average.coverage) # Median coverage of points
      medcount[i] <- length(medgc$average.coverage)
    }
    # Perform least squares regression
    reg <- lm(medy[medi] ~ poly(medx[medi], 3, raw=TRUE), weights=medcount[medi])
    p0 <- summary(reg)$coefficients[1,1]
    p1 <- summary(reg)$coefficients[2,1]
    p2 <- summary(reg)$coefficients[3,1]
    p3 <- summary(reg)$coefficients[4,1]
  }
  if(itercount>=100)
    print("Parameters to regression polynomial did not converge. Proceeding with final parameter values...")
  #print(paste("regression polynomial from inferred ploidy", reg))

  coverage
}

.normalizeRunningMedian <- function(coverage) {
  # Calculate here the average coverage across all regions
  totlen <- sum(coverage$length,na.rm=TRUE)
  totcov <- sum(coverage$length*coverage$average.coverage,na.rm=TRUE) # Total coverage across all regions
  avgcov <- totcov/totlen # Average coverage - used post-normalization

  # Begin running median normalization here
  # Get range of GC values to normalize over
  ############# medKP replace with coverage$diploid
  pmin <- as.integer(min(coverage[which(coverage$usable==TRUE),]$gc_bias)*1000)/1000 - 0.005
  pmax <- as.integer(max(coverage[which(coverage$usable==TRUE),]$gc_bias)*1000)/1000 + 0.005
  imin <- as.integer(min(coverage[which(coverage$diploid==TRUE),]$gc_bias)*1000)/1000 - 0.005
  imax <- as.integer(max(coverage[which(coverage$diploid==TRUE),]$gc_bias)*1000)/1000 + 0.005

  #Normalize
  coverage$NormalizedCoverage <- 0
  medDiploid <- data.frame()
  for(i in seq(pmin,pmax,0.01)) {
    # Big window limits
    i_low  <- i - 0.05
    i_high <- i + 0.05
    # Smaller window limits
    interval <- (i_high-i_low)/10.0
    medvals <- c(rep(NA,10))
    ct <- 1

    for(j in seq(i_low, i_high, interval)) { # Compute median of all small windows
          if(i>=imin & i<=imax)
        medvals[ct] <- median(coverage$average.coverage[which(coverage$diploid==TRUE & j<=coverage$gc_bias & coverage$gc_bias<j+0.01)], na.rm=TRUE)
          else
        medvals[ct] <- median(coverage$average.coverage[which(coverage$usable==TRUE & j<=coverage$gc_bias & coverage$gc_bias<j+0.01)], na.rm=TRUE)
      ct <- ct + 1
    }
    denom <- median(medvals, na.rm=TRUE) # Median of small window medians

    # Assign normalized coverage as: (original coverage/running median)*Average coverage across all regions
    coverage[which(i-0.005<=coverage$gc_bias & coverage$gc_bias<i+0.005 & coverage$usable==TRUE),"NormalizedCoverage"] <- avgcov*coverage[which(i-0.005<=coverage$gc_bias & coverage$gc_bias<i+0.005 & coverage$usable==TRUE),"average.coverage"] / denom
    gcNum <- median(coverage[which(i-0.005<=coverage$gc_bias & coverage$gc_bias<i+0.005 & coverage$usable==TRUE),"NormalizedCoverage"], na.rm=TRUE)
    imedDiploid <- cbind(i,denom,gcNum) # Pre and post-normalization medians
    medDiploid <- rbind(medDiploid,imedDiploid)
  }
  coverage[which(coverage$NormalizedCoverage==0),"NormalizedCoverage"] <- coverage[which(coverage$NormalizedCoverage==0),"average.coverage"]

  coverage$probe <- coverage$Target
  coverage$average.coverage <- coverage$NormalizedCoverage
  coverage$coverage <- coverage$average.coverage * coverage$length
  coverage$Target <- NULL
  colnames(medDiploid)[1] <- "gcIndex"

  ret <- list(coverage = coverage, medDiploid = medDiploid)

  ret
}

