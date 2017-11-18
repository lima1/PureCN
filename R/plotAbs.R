#' Plots for analyzing PureCN solutions
#' 
#' This function provides various plots for finding correct purity and ploidy
#' combinations in the results of a \code{\link{runAbsoluteCN}} call.
#' 
#' 
#' @param res Return object of the \code{\link{runAbsoluteCN}} function.
#' @param ids Candidate solutions to be plotted. \code{ids=1} will draw the
#' plot for the maximum likelihood solution.
#' @param type Different types of plots. \code{hist} will plot a histogram,
#' assigning log-ratio peaks to integer values. \code{overview} will plot all
#' local optima, sorted by likelihood. \code{BAF} plots
#' something like a B-allele frequency plot known from SNP arrays: it plots
#' allele frequencies of germline variants (or most likely germline when status
#' is not available) against copy number. \code{AF} plots observed allelic
#' fractions against expected (purity), maximum likelihood (optimal
#' multiplicity) allelic fractions. 
#' \code{all} plots all, and is useful
#' for generate a PDF for a sample for manual inspection. 
#' @param chr If \code{NULL}, show all chromosomes, otherwise only the ones
#' specified (\code{type="BAF"} only).
#' @param germline.only If \code{TRUE}, show only variants most likely being
#' germline in BAF plot. Useful to set to \code{FALSE} (in combination with
#' \code{chr}) to study potential artifacts.
#' @param show.contour For \code{type="overview"}, display contour plot.
#' @param purity Display expected integer copy numbers for purity, defaults to
#' purity of the solution (\code{type="hist"} only).
#' @param ploidy Display expected integer copy numbers for ploidy, defaults to
#' ploidy of the solution (\code{type="hist"} only).
#' @param alpha Add transparency to the plot if VCF contains many variants
#' (>2000, \code{type="AF"} and \code{type="BAF"} only).
#' @param show.segment.means Show segment means in germline allele frequency
#' plot?  If \code{both}, show SNVs and segment means. If \code{SNV} show all
#' SNVs. Only for \code{type="AF"}.
#' @param max.mapping.bias Exclude variants with high mapping bias from
#' plotting. Note that bias is reported on an inverse scale; a variant with
#' mapping bias of 1 has no bias. (\code{type="AF"} and \code{type="BAF"}
#' only).
#' @param palette.name The default \code{RColorBrewer} palette.
#' @param \dots Additonal parameters passed to the \code{plot} function.
#' @return Returns \code{NULL}.
#' @author Markus Riester
#' @seealso \code{\link{runAbsoluteCN}}
#' @examples
#' 
#' data(purecn.example.output)
#' plotAbs(purecn.example.output, type="overview")
#' # plot details for the maximum likelihood solution (rank 1)
#' plotAbs(purecn.example.output, 1, type="hist")
#' plotAbs(purecn.example.output, 1, type="BAF")
#' plotAbs(purecn.example.output, 1, type = "BAF", chr="chr2")
#' 
#' @export plotAbs
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette adjustcolor
#' @importFrom graphics abline axis boxplot contour hist image
#'             legend lines par plot text mtext polygon points
#'             rect strwidth symbols barplot
#' @importFrom ggplot2 geom_boxplot geom_hline labs
plotAbs <- function(res, ids = NULL, 
type = c("hist", "overview", "BAF", "AF", "all"),
chr = NULL, germline.only = TRUE, show.contour = FALSE, purity = NULL, 
ploidy = NULL, alpha = TRUE, show.segment.means = c("SNV", "segments", "both"),
max.mapping.bias = 0.8, palette.name = "Paired", ...) {
    chr.hash <- res$input$chr.hash
    if (is.null(chr.hash)) {
        chr.hash <- .getChrHash(gsub(":.*$","",res$input$log.ratio[,1]))
    }    
    centromeres <- .getCentromeres(res)

    sd.seg <- ifelse(is.null(res$input$log.ratio.sdev), 0.4, 
        res$input$log.ratio.sdev)

    .ar <- function(C) {
        (purity * C + 2*(1-purity))/( purity*ploidy + 2*(1-purity))  
    }
    .ardnorm <- function(x,C) dnorm(x, mean=log2(.ar(C)), sd=sd.seg, log=TRUE)
    .idealPeaks <- function(seg) {
        seg.split <- split(seg, seg$C)
        idx <- sapply(seg.split, nrow) > 2 | 
            ( sapply(seg.split, function(x) x$C[1]) <= 7 &  
              sapply(seg.split, function(x) x$C[1]) %in% 0:7 )
        peak.ideal.means <- log2(sapply(sapply(seg.split[idx], 
            function(x) x$C[1]), function(j) .ar(j)))
    }

    type <- match.arg(type)
    if (type=="hist") {
        if (is.null(ids)) ids <- seq_along(res$results)
        for (i in ids) {
            par.mar <- par("mar")
            par(mar=c(c(5, 4, 5, 2) + 0.1))
            seg <- res$results[[i]]$seg

            custom.solution <- TRUE
            if (is.null(purity)) { 
                purity <- res$results[[i]]$purity
            }    
            if (is.null(ploidy)) {
                ploidy <- res$results[[i]]$ploidy
                custom.solution <- FALSE
            }    
            peak.ideal.means <- .idealPeaks(seg)
            if (custom.solution) {
                peak.ideal.means <- log2(sapply(0:7, function(j) .ar(j)))
                names(peak.ideal.means) <- 0:7
            }    

            main <- paste("Purity:", round(purity, digits=2), 
                          " Tumor ploidy:", round( ploidy, digits=3))

            if (!is.null(res$results[[i]]$bootstrap.value)) {
                main <- paste(main, " Bootstrap value:", 
                    round(res$results[[i]]$bootstrap.value, digits=2))
            }    
            
            # try to avoid plotting log-ratio outliers from artifacts
            minLogRatio <- max(min(seg$seg.mean), peak.ideal.means[1]*2)
            logRatio <- do.call(c, lapply(seq_len(nrow(seg)), function(i)
                                rep(seg$seg.mean[i], seg$num.mark[i])))
            logRatio <- logRatio[logRatio>minLogRatio]

            h <- hist(logRatio,breaks=75, plot=FALSE)
            h$density <- h$counts/sum(h$counts)
            plot(h, freq=FALSE, xlab="log2 ratio",
                ylab="Fraction Genome", main=main, col="#377EB8", ...)
            abline(v=peak.ideal.means, lty=3, lwd=2, col="darkgrey")
            axis(side=3, at=peak.ideal.means, 
                labels=names(peak.ideal.means), tick=FALSE, padj=1)
            par(par.mar)

#zz<- lapply(0:7, function(C) curve(.ardnorm(x,C),log2(.ar(C))-0.5, log2(.ar(C))+0.5, xlim=c(-1,1.5), add=C>0))

        }
    } else if (type=="BAF") {
        if (is.null(res$input$vcf)) {
            .stopUserError("runAbsoluteCN was run without a VCF file.")
        }
        if (is.null(ids)) ids <- seq_along(res$results)
        for (i in ids) {
            r <- .getVariantPosteriors(res, i, NULL)
            if (is.null(r)) next
            seg <- res$results[[i]]$seg

            purity <- res$results[[i]]$purity
            ploidy <- res$results[[i]]$ploidy
            peak.ideal.means <- .idealPeaks(seg)
            GoF <- res$results[[i]]$GOF
            if (is.null(GoF)) {
                GoF <- .getGoF(res$results[[i]])
            }    
            r2 <- paste0(round(GoF * 100, digits=1),"%")

            par(mfrow=c(3,1))
            par(mar=c(c(5, 4, 5, 2) + 0.1))
            main <- paste(
                "Purity:", round(res$results[[i]]$purity[[1]], digits=2), 
                " Tumor ploidy:", round( res$results[[i]]$ploidy, digits=3), 
                " SNV log-likelihood:", 
                    round(res$results[[i]]$SNV.posterior$llik, digits=2),
                " GoF:", r2,    
                " Mean coverage:", 
                    paste(round(apply(geno(res$input$vcf)$DP,2,mean)),
                    collapse=";") )

            if (is.null(chr)) {
                if (germline.only) {
                    idx <-!r$ML.SOMATIC
                } else {
                    idx <- rep(TRUE, nrow(r))
                } 
            } else {
                idx <-r$chr %in% chr
                if (!sum(idx)) {
                    .stopUserError(paste(chr, collapse=","), 
                    " not valid chromosome name(s). ",
                    "Valid names are: ",  
                    paste(unique(r$chr), collapse=","))
                }
                if (germline.only) idx[r$ML.SOMATIC[idx]] <- FALSE
            }    
            if (!is.null(max.mapping.bias)) {
                idx <- idx & r$MAPPING.BIAS >= max.mapping.bias
            }
            #mycol <- ifelse(.strip.chr.name(r[idx,1], chr.hash) %% 2, 
            #    "#E41A1C", "#377EB8")
            mycol <- "#377EB8"
            mypch <- ifelse(r$GERMLINE.CONTHIGH > 0.5, 2, 
                ifelse(r$GERMLINE.CONTLOW>0.5, 3, 1))[idx]
            myalpha <- ifelse(alpha && nrow(r) > 2000, 2000/nrow(r), 1)

            tmp <- data.frame(chr=levels(r[,1]), start=match(levels(r[,1]), r[idx, 1]), 
                stringsAsFactors=FALSE)
            tmp <- tmp[complete.cases(tmp),,drop=FALSE]
            tmp <- tmp[order(.strip.chr.name(tmp$chr,chr.hash)),]
            tmp$end <- c(tmp[-1,2]-1, nrow(r))
            cids <- NULL
            if (!is.null(centromeres)) {
                cids <- sapply(tmp$chr, function(x) { 
                    id <- start(centromeres)[which(seqnames(centromeres)==x)]; 
                    which(r$chr[idx]==x & r$start[idx] >= id)[1] })
            }
            segment.log.ratio <- res$results[[i]]$seg$seg.mean[r$seg.id]
            segment.log.ratio.lines <- .toLines(ss=segment.log.ratio[idx])
            segment.M.log.ratio.lines <- .toLines(
                peak.ideal.means[as.character(
                res$results[[i]]$SNV.posterior$posteriors$ML.M.SEGMENT[idx]
            )])
            
            # calculate expected segment B-allelic fractions
            purity <- res$results[[i]]$purity
            ploidy <- res$results[[i]]$ploidy
            b1 <- ((purity*r$ML.M.SEGMENT[idx])+(1-purity))/
                ((purity*r$ML.C[idx])+2*(1-purity))
            b2 <- 1-b1
            segment.b1.lines <- .toLines(ss=b1)
            segment.b2.lines <- .toLines(ss=b2)

            if (!is.null(chr) && length(chr)==1) {
                x <- r$start[idx]/1000
                logRatio <- res$input$log.ratio
                logRatio <- logRatio[seqnames(logRatio) %in% chr]
                xLogRatio <- start(logRatio)/1000
                plot(x, r$AR[idx],ylab="B-Allele Frequency", 
                    xlab="Pos (kbp)",main=paste(main, " Chromosome:", chr), 
                    col=mycol, pch=mypch, xlim=range(xLogRatio), ...)
                centromerePos <- NULL
                if (!is.null(res$input$centromere)) {
                    centromerePos <- start(centromeres)[
                        which(seqnames(centromeres)==chr)]/1000
                }    

                segment.b1.lines[,1] <- x[segment.b1.lines[,1]]
                segment.b2.lines[,1] <- x[segment.b2.lines[,1]]
                lines(segment.b1.lines, col="black", lwd=3)
                lines(segment.b2.lines, col="black", lwd=3)

                abline(h=0.5, lty=3, col="grey")
                abline(v=centromerePos, lty=3, col="grey")
                 
                plot(xLogRatio, logRatio$log.ratio, ylab="Copy Number log-ratio", 
                    xlab="Pos (kbp)", 
                    col=adjustcolor("grey", alpha.f=ifelse(myalpha<1,0.75,1)), 
                    ylim=range(segment.log.ratio*1.1),
                    ... )
                points(x, r$log.ratio[idx], col=mycol, pch=mypch)
                segment.log.ratio.lines[,1] <- x[segment.log.ratio.lines[,1]]
                lines(segment.log.ratio.lines, col="black", lwd=3)
                segment.M.log.ratio.lines[,1] <- x[segment.M.log.ratio.lines[,1]]
                lines(segment.M.log.ratio.lines, col="grey", lwd=3)
                abline(h=0, lty=3, col="grey")
                abline(h=peak.ideal.means, lty=2, col="grey")
                axis(side=4,at=peak.ideal.means, 
                    labels=names(peak.ideal.means))
                abline(v=centromerePos, lty=3, col="grey")
                plot(x, r$ML.C[idx], ylab="Maximum Likelihood Copy Number", 
                    xlab="Pos (kbp)", xlim=range(xLogRatio),
                    ylim=c(0,min(7, max(r$ML.C[!r$ML.SOMATIC]))), ... )
                points(x, r$ML.M.SEGMENT[idx], col="grey")
                abline(v=centromerePos, lty=3, col="grey")
            } else {
                plot(r$AR[idx],ylab="B-Allele Frequency", xlab="SNV Index",
                    main=main, type="n", ...)
                rect(tmp$start, par("usr")[3], tmp$end+1, par("usr")[4], 
                     col=ifelse(seq(nrow(tmp))%%2, "#deebf7", "white"), border=NA)
                rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4])
                points(r$AR[idx], col=adjustcolor(mycol, alpha.f=myalpha), pch=mypch)
                lines(segment.b1.lines, col="black", lwd=3)
                lines(segment.b2.lines, col="black", lwd=3)
                axis(side=3, at=(tmp[,3]+tmp[,2])/2, 
                    labels=.strip.chr.name(tmp[,1], chr.hash), 
                    tick=FALSE, padj=1)
                abline(h=0.5, lty=3, col="grey")
                #abline(v=tmp[,2], lty=3, col="grey")
                abline(v=cids, lty=3, col="grey")
                main <- paste("SCNA-fit Log-Likelihood:", 
                    round(res$results[[i]]$log.likelihood, digits=2) )

                myylim <- quantile(subset(r$log.ratio,
                    !is.infinite(r$log.ratio)), p=c(0.001, 1-0.001),na.rm=TRUE)
                myylim[1] <- floor(myylim[1])
                myylim[2] <- ceiling(myylim[2])

                plot(r$log.ratio[idx], ylab="Copy Number log-ratio", 
                    xlab="SNV Index", 
                    main=main, ylim=myylim, type="n",... )
                rect(tmp$start, par("usr")[3], tmp$end+1, par("usr")[4], 
                     col=ifelse(seq(nrow(tmp))%%2, "#deebf7", "white"), border=NA)
                rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4])
                points(r$log.ratio[idx], col=adjustcolor(mycol, alpha.f=myalpha), pch=mypch)
                lines(segment.log.ratio.lines, col="black", lwd=3)
                lines(segment.M.log.ratio.lines, col="grey", lwd=3)

                abline(h=0, lty=3, col="grey")
                abline(v=cids, lty=3, col="grey")
                abline(h=peak.ideal.means, lty=2, col="grey")
                axis(side=4,at=peak.ideal.means, 
                    labels=names(peak.ideal.means))

                plot(r$ML.M.SEGMENT[idx], ylab="Maximum Likelihood Copy Number", 
                    xlab="SNV Index", 
                    ylim=c(0,min(7, max(r$ML.C[!r$ML.SOMATIC]))), col="grey",
                    type="n", ... )
                rect(tmp$start, par("usr")[3], tmp$end+1, par("usr")[4], 
                     col=ifelse(seq(nrow(tmp))%%2, "#deebf7", "white"), border=NA)
                rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4])
                points(r$ML.M.SEGMENT[idx], col="grey")
                points(r$ML.C[idx], col="black")
                abline(v=cids, lty=3, col="grey")
            } 
        }
    } else if (type=="AF") {
        if (is.null(res$input$vcf)) {
            .stopUserError("runAbsoluteCN was run without a VCF file.")
        }
        if (is.null(ids)) ids <- seq_along(res$results)
        for (i in ids) {
            r <- .getVariantPosteriors(res, i, max.mapping.bias)
            if (is.null(r)) next

            vcf <- res$input$vcf[res$results[[i]]$SNV.posterior$vcf.ids]
            # brwer.pal requires at least 3 levels

            r <- .getAFPlotGroups(r, is.null(info(vcf)$SOMATIC))
            r <- r[order(r$group),]

            tmp <- I(brewer.pal(max(nlevels(r$group),3),
                name=palette.name ))

            mycol.palette <- data.frame(
                group=levels(r$group), 
                color=tmp[seq_len(nlevels(r$group))] 
            )

            mycol.palette$pch <- seq_len(nrow(mycol.palette))

            mycol <- mycol.palette$color[
                match(as.character(r$group), mycol.palette$group)]
            mypch <- mycol.palette$pch[
                match(as.character(r$group), mycol.palette$group)]

            main.color <- sort(table(mycol), decreasing=TRUE)[1]

            par(mfrow=c(2,2))
            main <- paste(
                "Purity:", round(res$results[[i]]$purity[[1]], digits=2), 
                " Tumor ploidy:", round( res$results[[i]]$ploidy, digits=3)
            )
            
            myalpha <- ifelse(alpha && nrow(r) > 2000, 2000/nrow(r), 1)
            myalpha <- 1
            mycex <- log10(r$depth*r$MAPPING.BIAS)-1

            seg <- res$results[[i]]$seg
            if (is.null(purity)) { 
                purity <- res$results[[i]]$purity
            }    
            if (is.null(ploidy)) {
                ploidy <- res$results[[i]]$ploidy
                custom.solution <- FALSE
            }    

            mylogratio.xlim <- quantile(subset(r$log.ratio,
                    !is.infinite(r$log.ratio)), p=c(0.001, 1-0.001),na.rm=TRUE)

            peak.ideal.means <- .idealPeaks(seg)
            scatter.labels <- paste0(r$ML.C,"m", r$ML.M)[!r$ML.SOMATIC]
            idx.labels <- !duplicated(scatter.labels) & 
                as.character(r$ML.C[!r$ML.SOMATIC]) %in% 
                names(peak.ideal.means)
            segment.means <- match.arg(show.segment.means)

            idx.nna <- !is.na(r$ML.M ) & !r$ML.SOMATIC
            y <- sapply(split(r$AR[idx.nna],  paste(r$ML.M, r$seg.id)[idx.nna]), mean)
            x <- sapply(split(r$log.ratio[idx.nna],  paste(r$ML.M, r$seg.id)[idx.nna]), mean)
            mlm <- sapply(split(r$ML.M[idx.nna],  paste(r$ML.M, r$seg.id)[idx.nna]), mean)
            mlc <- sapply(split(r$ML.C[idx.nna],  paste(r$ML.M, r$seg.id)[idx.nna]), mean)
            size <- sapply(split(r$ML.M[idx.nna],  paste(r$ML.M, r$seg.id)[idx.nna]), length)
            
            if (segment.means %in% c("both", "SNV")) {    
                plot(r$log.ratio[!r$ML.SOMATIC], r$AR[!r$ML.SOMATIC], 
                    col=adjustcolor(mycol[!r$ML.SOMATIC], alpha.f=myalpha), pch=mypch[!r$ML.SOMATIC], 
                    xlab="Copy Number log-ratio", ylab="Allelic fraction (germline)",
                    cex=mycex[!r$ML.SOMATIC],
                    xlim=mylogratio.xlim
                    )
                if (segment.means == "both") points(x,y,
                   col=adjustcolor("orange", alpha.f=0.2),pch=20,cex=log2(size))
            } else {
                plot(x,y,
                    col=adjustcolor("orange", alpha.f=0.5),pch=20,cex=log2(size),
                    xlab="Copy Number log-ratio", ylab="Allelic fraction (germline)",
                    xlim=mylogratio.xlim
                    )
            } 
            text(
                x=peak.ideal.means[
                    as.character(r$ML.C[!r$ML.SOMATIC])][idx.labels], 
                y=r$ML.AR[!r$ML.SOMATIC][idx.labels], 
                labels=scatter.labels[idx.labels]
            )
            plot(r$AR, log2(r$depth),col=adjustcolor(mycol, alpha.f=myalpha),
                pch=mypch, cex=mycex,
                xlab="Allelic fraction", ylab="Coverage (log2)")

            if (sum(r$ML.SOMATIC)>0) {

                scatter.labels <- paste0(r$ML.C,"m", r$ML.M)[r$ML.SOMATIC]

                idx.labels <- !duplicated(scatter.labels) & 
                    as.character(r$ML.C[r$ML.SOMATIC]) %in% 
                    names(peak.ideal.means)

                plot(r$log.ratio[r$ML.SOMATIC], r$AR[r$ML.SOMATIC], 
                    col=mycol[r$ML.SOMATIC], pch=mypch[r$ML.SOMATIC], 
                    cex=mycex[r$ML.SOMATIC],
                    xlab="Copy Number log-ratio", 
                    ylab="Allelic fraction (somatic)",
                    xlim=mylogratio.xlim
                    )
                legend("topright", legend=as.character(mycol.palette$group),
                    col=mycol.palette$color, 
                    pch=mycol.palette$pch, cex=0.8)
                
                estimatedContRate <- 
                    res$results[[i]]$SNV.posterior$posterior.contamination
                if (!is.null(estimatedContRate) && 
                    estimatedContRate > min(r$AR[r$ML.SOMATIC])) {
                    abline(h=estimatedContRate, col="red")
                    text(x=mylogratio.xlim[1], y=estimatedContRate+0.02, 
                        labels="Contamination", col="red", pos=4)
                }
                text(x=peak.ideal.means[
                    as.character(r$ML.C[r$ML.SOMATIC])][idx.labels], 
                    y=r$ML.AR[r$ML.SOMATIC][idx.labels], 
                    labels=scatter.labels[idx.labels])
                idxSomatic <- !grepl("germline|dbSNP|contamination", as.character(r$group))
                if (sum(idxSomatic)) {
                    colSomatic <- mycol.palette$color[match(names(sort(table(r$group[idxSomatic]), 
                        decreasing=TRUE)[1]), mycol.palette$group)]
                    hist(r$CELLFRACTION[idxSomatic], col=colSomatic,
                        xlab="Cellular fraction", main="")
                    cellFractions <- seq(0.01,1,0.01)
                    power <- sapply(cellFractions, function(pw) 
                        calculatePowerDetectSomatic(mean(r$depth, na.rm=TRUE), 
                            purity=res$results[[i]]$purity,
                            ploidy=res$results[[i]]$ploidy, cell.fraction=pw,
                            verbose=FALSE)$power)
                    if (max(power) > 0.8) {
                        minFraction <- cellFractions[which(power > 0.8)[1]]
                        abline(v=minFraction)
                        axis(side=3, at=minFraction, 
                            labels="Power 0.8", tick=FALSE, padj=1)
                    } else {
                        abline(v=1)
                        axis(side=3, at=1, 
                            labels=paste("Power", round(max(power),digits=2)), 
                            tick=FALSE, padj=1)
                    }    
                }
            } else {
                legend("bottomright", legend=as.character(mycol.palette$group),
                    col=mycol.palette$color, pch=mycol.palette$pch)
            }
        }
    } else if (type == "all") {
        if (is.null(ids)) ids <- seq_along(res$results)
        for (i in ids) {
            par(mfrow=c(1,1))
            plotAbs(res,i, type="hist")
            if (!is.null(res$input$vcf)) {
                plotAbs(res,i, type="BAF",...)
                plotAbs(res,i, type="AF",...)
            }
        }    
    } else {
        mycol <-  ifelse(sapply(res$results, function(x) x$flag), "yellow", 
            "white")
        mycolBg <-  ifelse(sapply(res$results, function(x) x$flag), "black", 
            "black")
        myfont <-  ifelse(sapply(res$results, function(x) x$flag), 1, 1)
        main <- NULL
        parm <- par("mar")
        par(mar= c(5, 4, 4, 4) + 0.1)
        #if ("darkgrey" %in% myfont) main="Italics: SCNA-fitting flagged."
        xc <- .matrixTotalPloidyToTumorPloidy(res$candidates$all)
        xc[is.infinite(xc)] <- min(xc[!is.infinite(xc)])
        xc[xc < quantile(xc, p=0.2)] <- quantile(xc, p=0.2)
        
        mycol.image <- colorRampPalette(rev(brewer.pal(n = 7, 
            name = "RdYlBu")))(100)
        image(as.numeric(colnames(xc)), as.numeric(rownames(xc)), 
            t(xc)-max(xc), col=mycol.image, xlab = "Purity", 
            ylab = "Ploidy",main = main,...)
        .legend.col(col=mycol.image, lev=min(xc):max(xc), 
            ylim=quantile(as.numeric(rownames(xc)), p=c(0,1)))

        if (show.contour) contour(as.numeric(colnames(xc)), 
            as.numeric(rownames(xc)), t(xc), add=TRUE)

        symbols(x=sapply(res$results, function(x) min(0.95, x$purity)) - 0.02, 
            y=sapply(res$results, function(x) x$ploidy) - 0.1, 
            circles=rep(mean(sapply(seq_along(res$results), strwidth, 
                cex=1.25)), length(res$results)), 
            bg=mycolBg, inches=FALSE, add=TRUE)

        text( sapply(res$results, function(x) min(0.95, x$purity))-0.02,
            sapply(res$results, function(x) x$ploidy)-0.1, 
            seq_along(res$results), col=mycol, cex=1.2,font=myfont)
        par(mar=parm)
    }
}

.toLines <- function(
### "segments" already segmented log-ratios into a list for plotting
ss) {
    #no lines to draw
    if (length(ss)<2) return(matrix(1))
    # avoid the corner case when last SNV is in different segment
    if (length(ss)>2) ss[length(ss)] <- ss[length(ss)-1]

    #get breakpoints
    bp <- sapply(2:(length(ss)),function(i) ss[i] != ss[i-1])
    bp <- c(TRUE, bp)
    xy.start <- cbind(x=which(bp), y=ss[which(bp)])
    xy.end <- cbind(x=c(x=which(bp)[-1], length(bp)), y=ss[which(bp)])
    xy.end.na <- cbind(x=c(x=which(bp)[-1], length(bp)), y=NA)
    xy <- rbind(xy.end, xy.end.na, xy.start)
    xy <- xy[order(xy[,1]),]
}

.matrixTotalPloidyToTumorPloidy <- function(ca) {
    ploidy <- as.numeric(rownames(ca))

    tumor.ploidy <- lapply(as.numeric(colnames(ca)), function(purity) 
        (ploidy - 2*(1-purity))/purity)

    tumor.ploidy.i <- lapply(tumor.ploidy, function(x) 
        sapply(ploidy, function(i) ifelse(min(abs(x-i)) < 1,
             which.min(abs(x-i) ), NA)))

    cc <- do.call(cbind, lapply(seq_along(tumor.ploidy.i), function(i) 
        ca[tumor.ploidy.i[[i]],i]))

    colnames(cc) <- colnames(ca)
    rownames(cc) <- rownames(ca)
    cc
}

.legend.col <- function(col, lev, ylim){
    n <- length(col)
    bx <- par("usr")
     
    box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
    bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
    box.cy <- c(bx[3], bx[3])
    box.sy <- (bx[4] - bx[3]) / n
     
    xx <- rep(box.cx, each = 2)

    pxpd = par("xpd")
    par(xpd = TRUE)
    for(i in seq_len(n)){
        yy <- c(box.cy[1] + (box.sy * (i - 1)),
        box.cy[1] + (box.sy * (i)),
        box.cy[1] + (box.sy * (i)),
        box.cy[1] + (box.sy * (i - 1)))
        polygon(xx, yy, col = col[i], border = col[i])
    }
    axis(side=4, at=c(ylim[1], (ylim[2]+ylim[1])/2, ylim[2]), tick=FALSE, 
    labels=round(c(min(lev), median(lev), max(lev))))
    mtext("Copy number Log-Likelihood",side=4,line=2)
    par(xpd = pxpd)
}

.getVariantPosteriors <- function(res, i, max.mapping.bias=NULL) {
    r <- res$results[[i]]$SNV.posterior$posteriors
    if (!is.null(r) && !is.null(max.mapping.bias)) {
        r <- r[r$MAPPING.BIAS >= max.mapping.bias,]
    }    
    r
}
   
.getAFPlotGroups <- function(r, single.mode) {
    if (single.mode) {
        groupLevels <- c("dbSNP/germline", "dbSNP/somatic", "novel/somatic",
            "novel/germline", "COSMIC/germline", "COSMIC/somatic", "contamination") 
        r$group <- groupLevels[1]
        r$group[r$prior.somatic < 0.1 & r$ML.SOMATIC] <- groupLevels[2]
        r$group[r$prior.somatic >= 0.1 & r$ML.SOMATIC] <- groupLevels[3]
        r$group[r$prior.somatic >= 0.1 & !r$ML.SOMATIC] <- groupLevels[4]
        r$group[r$prior.somatic >= 0.9 & !r$ML.SOMATIC] <- groupLevels[5]
        r$group[r$prior.somatic >= 0.9 & r$ML.SOMATIC] <- groupLevels[6]
        r$group[r$GERMLINE.CONTLOW > 0.9 | 
            r$GERMLINE.CONTHIGH > 0.9] <- groupLevels[7]
        r$group <- factor(r$group, levels=groupLevels)
        return(r)
    } 
    groupLevels <- c("germline", "germline/ML somatic", "somatic", 
        "somatic/ML germline", "contamination")
    r$group <- groupLevels[1]
    r$group[r$prior.somatic < 0.1 & r$ML.SOMATIC] <- groupLevels[2]
    r$group[r$prior.somatic >= 0.1 & r$ML.SOMATIC] <- groupLevels[3]
    r$group[r$prior.somatic >= 0.1 & !r$ML.SOMATIC] <- groupLevels[4]
    r$group[r$GERMLINE.CONTLOW > 0.9 | 
        r$GERMLINE.CONTHIGH > 0.9] <- groupLevels[5]
    r$group <- factor(r$group, levels=groupLevels)
    r
}            

.plotLogRatios <- function(log.ratio, on.target) {
    containsOfftarget <- sum(on.target)!=length(on.target)
    if (!containsOfftarget) return(NULL)
    myylim <- quantile(subset(log.ratio,
        !is.infinite(log.ratio)), p=c(0.0001, 1-0.0001),na.rm=TRUE)
    plot(log.ratio, col=ifelse(on.target, "black", "red"), 
        pch=".",cex=3, ylim=myylim, ylab="log2 ratio")
    legend("bottomleft", legend=c("On-Target", "Off-Target"), ncol=2, fill=c("black", "red"))
}
    
