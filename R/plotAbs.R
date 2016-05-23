plotAbs <-
structure(function(# Plots for analyzing PureCN solutions
### This function provides various plots for finding correct 
### purity and ploidy combinations in the results of a runAbsoluteCN call.
res, 
### Return object of the runAbsoluteCN() function.
ids=NULL, 
### Candidate solutions to be plotted. ids=1 will draw the 
### plot for the maximum likelihood solution.
type=c("hist", "overview", "BAF", "AF", "LOH", "all"),
### Different types of plots. "hist" will plot a histogram, 
### assigning log-ratio peaks to integer values. "overview" will plot all 
### local optima, sorted by likelihood. "BAF" plots something like a B-allele
### frequency plot known from SNP arrays: it plots allele frequencies of 
### germline variants (or most likely germline when status is not available) 
### against copy number. AF plots observed allelic fractions against expected
### (purity), maximum likelihood (optimal multiplicity) allelic fractions. 
### "all" plots all, and is useful for generate a PDF for a sample for manual
### inspection.
chr=NULL,
### If NULL, show all chromosomes, otherwise only the ones 
### specified (type=BAF only).
germline.only=TRUE,
### If TRUE, show only variants most likely being germline in 
### BAF plot. Useful to set to FALSE (in combination with chr) to study 
### potential artifacts.
show.contour=FALSE,
### For type overview, display contour plot.
purity=NULL,
### Display expected integer copy numbers for purity, defaults 
### to purity of the solution (type=hist only).
ploidy=NULL,
### Display expected integer copy numbers for ploidy, defaults 
### to ploidy of the solution (type=hist only).
alpha=TRUE,
### Add transparency to the plot if VCF contains many variants 
### (>1000, type=AF only). 
...
### Additonal parameters passed to the plot() function. 
) {
    sd.seg <- ifelse(is.null(res$input$log.ratio.sdev), 0.4, res$input$log.ratio.sdev)

    .ar <- function(C) { 
        (purity * C + 2*(1-purity))/( purity*ploidy + 2*(1-purity))  
    }
    .ardnorm <- function(x,C) dnorm(x, mean=log2(.ar(C)), sd=sd.seg, log=TRUE)

    type <- match.arg(type)
    if (type=="hist") {
        if (is.null(ids)) ids <- 1:length(res$results)
        for (i in ids) {
            par.mar <- par("mar")
            par(mar=c(c(5, 4, 5, 2) + 0.1))
            seg <- res$results[[i]]$seg
            seg.split <- split(seg, seg$C)
            idx <- sapply(seg.split, nrow) > 2 | 
                ( sapply(seg.split, function(x) x$C[1]) <= 7 &  
                  sapply(seg.split, function(x) x$C[1]) %in% 0:7 )

            custom.solution <- TRUE
            if (is.null(purity)) { 
                purity <- res$results[[i]]$purity
            }    
            if (is.null(ploidy)) {
                ploidy <- res$results[[i]]$ploidy
                custom.solution <- FALSE
            }    

            peak.ideal.means <- log2(sapply(sapply(seg.split[idx], 
                function(x) x$C[1]), function(j) .ar(j)))
            if (custom.solution) {
                peak.ideal.means <- log2(sapply(0:7, function(j) .ar(j)))
                names(peak.ideal.means) <- 0:7
            }    

            main <- paste("Purity:", round(purity, digits=2), 
                          "Tumor Ploidy:", round( ploidy, digits=3))

            h <- hist(do.call(c, lapply(1:nrow(seg), function(i)
                    rep(seg$seg.mean[i], seg$num.mark[i]))),breaks=100, 
                    plot=FALSE)
            h$density <- h$counts/sum(h$counts)
            plot(h, freq=FALSE, xlab="log2 ratio",
                ylab="Fraction Genome", main=main,...)
            abline(v=peak.ideal.means, lty=3, lwd=2, col="darkgrey")
            axis(side=3, at=peak.ideal.means, 
                labels=names(peak.ideal.means), tick=FALSE, padj=1)
            par(par.mar)

#zz<- lapply(0:7, function(C) curve(.ardnorm(x,C),log2(.ar(C))-0.5, log2(.ar(C))+0.5, xlim=c(-1,1.5), add=C>0))

        }
    } else if (type=="BAF") {
        if (is.null(res$input$vcf)) stop("runAbsoluteCN was run without a VCF file.")
        if (is.null(ids)) ids <- 1:length(res$results)
        for (i in ids) {
            r <- res$results[[i]]$SNV.posterior$beta.model$posterior
            if (is.null(r)) next
            par(mfrow=c(3,1))
            par(mar=c(c(5, 4, 5, 2) + 0.1))
            main <- paste(
                "Purity:", round(res$results[[i]]$purity[[1]], digits=2), 
                " Tumor Ploidy:", round( res$results[[i]]$ploidy, digits=3), 
                " SNV Log-Likelihood:", 
                    round(res$results[[i]]$SNV.posterior$beta$llik, digits=2),
                " Mean Coverage:", 
                    paste(round(apply(geno(res$input$vcf)$DP,2,mean)),
                    collapse=";") )

            if (is.null(chr)) {
                if (germline.only) {
                    idx <-!r$ML.SOMATIC
                } else {
                    idx <- rep(TRUE, nrow(r))
                } 
            } else {
                idx <-r$seqnames %in% chr
                if (germline.only) idx[r$ML.SOMATIC[idx]] <- FALSE
            }    
            mycol <- ifelse(.strip.chr.name(r[idx,1]) %% 2, 
                "#E41A1C", "#377EB8")
            mypch <- ifelse(r$GERMLINE.CONTHIGH > 0.5, 2, 
                ifelse(r$GERMLINE.CONTLOW>0.5, 3, 1))[idx]

            tmp <- cbind(levels(r[,1]), match(levels(r[,1]), r[idx, 1]))
            tmp <- tmp[complete.cases(tmp),,drop=FALSE]

            segment.log.ratio <- res$results[[i]]$seg$seg.mean[
                res$results[[i]]$SNV.posterior$beta.model$segment.ids]
            segment.log.ratio.lines <- .toLines(ss=segment.log.ratio[idx])
            
            if (!is.null(chr) && length(chr)==1) {
                x <- r$start[idx]/1000
                plot(x, r$AR[idx],ylab="B-Allele Frequency", 
                    xlab="Pos (kbp)",main=paste(main, " Chromosome:", chr), 
                    col=mycol, pch=mypch, ...)
                abline(h=0.5, lty=3, col="grey")
                main <- paste("SCNA-fit Log-Likelihood:", 
                    round(res$results[[i]]$log.likelihood, digits=2) )
                plot(x, r$Log.Ratio[idx], ylab="Copy Number log-ratio", 
                    xlab="Pos (kbp)", col=mycol, main=main, pch=mypch,... )
                segment.log.ratio.lines[,1] <- x[segment.log.ratio.lines[,1]]
                lines(segment.log.ratio.lines, col="grey", lwd=3)
                abline(h=0, lty=3, col="grey")
                plot(x, r$ML.C[idx], ylab="Maximum Likelihood Copy Number", 
                    xlab="Pos (kbp)", 
                    ylim=c(0,min(7, max(r$ML.C[!r$ML.SOMATIC]))), ... )
            } else {
                plot(r$AR[idx],ylab="B-Allele Frequency", xlab="SNV Index",
                    main=main, col=mycol, pch=mypch, ...)
                axis(side=3, at=tmp[,2], labels=tmp[,1], tick=FALSE, padj=1)
                abline(h=0.5, lty=3, col="grey")
                abline(v=tmp[,2], lty=3, col="grey")
                main <- paste("SCNA-fit Log-Likelihood:", 
                    round(res$results[[i]]$log.likelihood, digits=2) )

                myylim <- quantile(subset(r$Log.Ratio,
                    !is.infinite(r$Log.Ratio)), p=c(0.001, 1-0.001),na.rm=TRUE)
                myylim[1] <- floor(myylim[1])
                myylim[2] <- ceiling(myylim[2])

                plot(r$Log.Ratio[idx], ylab="Copy Number log-ratio", 
                    xlab="SNV Index", col=mycol, main=main, pch=mypch, 
                    ylim=myylim, ... )
                lines(segment.log.ratio.lines, col="grey", lwd=3)

                abline(h=0, lty=3, col="grey")
                abline(v=tmp[,2], lty=3, col="grey")

                plot(r$ML.C[idx], ylab="Maximum Likelihood Copy Number", 
                    xlab="SNV Index", 
                    ylim=c(0,min(7, max(r$ML.C[!r$ML.SOMATIC]))), ... )
                abline(v=tmp[,2], lty=3, col="grey")
            } 
        }
    } else if (type=="AF") {
        if (is.null(res$input$vcf)) stop("runAbsoluteCN was run without a VCF file.")
        if (is.null(ids)) ids <- 1:length(res$results)
        for (i in ids) {
            r <- res$results[[i]]$SNV.posterior$beta$posterior
            if (is.null(r)) next
            vcf <- res$input$vcf[res$results[[i]]$SNV.posterior$beta$vcf.ids]
            # brwer.pal requires at least 3 levels
            tmp <- I(brewer.pal(max(nlevels(as.factor(r$Prior.Somatic)),3),
                name="Dark2" ))

            mycol.palette <- data.frame(
                priors=levels(as.factor(r$Prior.Somatic)), 
                color=tmp[1:nlevels(as.factor(r$Prior.Somatic))] 
            )

            mycol.palette$pch <- 1:nrow(mycol.palette)

            mycol <- mycol.palette$color[
                match(as.character(r$Prior.Somatic), mycol.palette$priors)]
            mypch <- mycol.palette$pch[
                match(as.character(r$Prior.Somatic), mycol.palette$priors)]

            main.color <- sort(table(mycol), decreasing=TRUE)[1]

            par(mfrow=c(2,2))
            main <- paste(
                "Purity:", round(res$results[[i]]$purity[[1]], digits=2), 
                " Tumor Ploidy:", round( res$results[[i]]$ploidy, digits=3)
            )
            if (!alpha || main.color < 1000) {
                plot(r$ML.AR[!r$ML.SOMATIC], r$AR[!r$ML.SOMATIC],  
                    col=mycol[!r$ML.SOMATIC], pch=mypch[!r$ML.SOMATIC], 
                    xlab="Expected allelic fraction", 
                    ylab="Allelic fraction (germline)", main=main,...)
            } else {
                smoothScatter(
                    r$ML.AR[!r$ML.SOMATIC & mycol==names(main.color)], 
                    r$AR[!r$ML.SOMATIC & !r$ML.SOMATIC & mycol==names(main.color)], 
                    colramp=colorRampPalette(c("white", names(main.color))),
                    xlab="Expected allelic fraction", 
                    ylab="Allelic fraction (germline)", main=main, 
                    transformation = function(x) x,...)
                points(r$ML.AR[!r$ML.SOMATIC & mycol != names(main.color)], 
                        r$AR[!r$ML.SOMATIC &  mycol != names(main.color)], 
                        col=mycol[!r$ML.SOMATIC &  mycol != names(main.color)], 
                        pch=mypch[!r$ML.SOMATIC &  mycol != names(main.color)])
            }    
            abline(a=0, b=1, lty=3, col="grey")

            seg <- res$results[[i]]$seg
            seg.split <- split(seg, seg$C)
            idx <- sapply(seg.split, nrow) > 2 | 
                ( sapply(seg.split, function(x) x$C[1]) <= 7 &  
                  sapply(seg.split, function(x) x$C[1]) %in% 0:7 )
            if (is.null(purity)) { 
                purity <- res$results[[i]]$purity
            }    
            if (is.null(ploidy)) {
                ploidy <- res$results[[i]]$ploidy
                custom.solution <- FALSE
            }    

            peak.ideal.means <- log2(sapply(sapply(seg.split[idx], 
                            function(x) x$C[1]), function(j) .ar(j)))
            scatter.labels <- paste(r$ML.C,"m", r$ML.M, sep="")[!r$ML.SOMATIC]
            idx.labels <- !duplicated(scatter.labels) & 
                as.character(r$ML.C[!r$ML.SOMATIC]) %in% 
                names(peak.ideal.means)
            if (!alpha || main.color < 1000) {
                plot(r$Log.Ratio[!r$ML.SOMATIC], r$AR[!r$ML.SOMATIC], 
                    col=mycol[!r$ML.SOMATIC], pch=mypch[!r$ML.SOMATIC], 
                    xlab="Copy Number log-ratio", ylab="Allelic fraction (germline)")
            } else {
                smoothScatter(
                    r$Log.Ratio[!r$ML.SOMATIC & mycol==names(main.color)], 
                    r$AR[!r$ML.SOMATIC & !r$ML.SOMATIC & mycol==names(main.color)], 
                    colramp=colorRampPalette(c("white", names(main.color))),
                    xlab="Copy Number log-ratio", ylab="Allelic fraction (germline)",
                    transformation = function(x) x
                    )
                points(r$Log.Ratio[!r$ML.SOMATIC & mycol != names(main.color)], 
                        r$AR[!r$ML.SOMATIC &  mycol != names(main.color)], 
                        col=mycol[!r$ML.SOMATIC &  mycol != names(main.color)], 
                        pch=mypch[!r$ML.SOMATIC &  mycol != names(main.color)])
            }    
            text(
                x=peak.ideal.means[
                    as.character(r$ML.C[!r$ML.SOMATIC])][idx.labels], 
                y=r$ML.AR[!r$ML.SOMATIC][idx.labels], 
                labels=scatter.labels[idx.labels]
            )

            if (sum(r$ML.SOMATIC)>0) {
                plot(r$ML.AR[r$ML.SOMATIC], r$AR[r$ML.SOMATIC], 
                    col=mycol[r$ML.SOMATIC], pch=mypch[r$ML.SOMATIC], 
                    xlab="Expected allelic fraction", 
                    ylab="Allelic fraction (somatic)",...)
                abline(a=0, b=1, lty=3, col="grey")
                legend("bottomright", legend=paste("Prior Somatic", 
                    round(as.numeric(as.character(mycol.palette$priors)),
                    digits=4)) , col=mycol.palette$color, 
                    pch=mycol.palette$pch, cex=0.8)

                scatter.labels <- paste(r$ML.C,"m", r$ML.M, 
                    sep="")[r$ML.SOMATIC]

                idx.labels <- !duplicated(scatter.labels) & 
                    as.character(r$ML.C[r$ML.SOMATIC]) %in% names(peak.ideal.means)

                plot(r$Log.Ratio[r$ML.SOMATIC], r$AR[r$ML.SOMATIC], 
                    col=mycol[r$ML.SOMATIC], pch=mypch[r$ML.SOMATIC], 
                    xlab="Copy Number log-ratio", 
                    ylab="Allelic fraction (somatic)")

                text(x=peak.ideal.means[
                    as.character(r$ML.C[r$ML.SOMATIC])][idx.labels], 
                    y=r$ML.AR[r$ML.SOMATIC][idx.labels], 
                    labels=scatter.labels[idx.labels])
            } else {
                legend("bottomright", legend=paste("Prior Somatic", 
                    round(as.numeric(as.character(mycol.palette$priors)),
                    digits=4)), col=mycol.palette$color, pch=mycol.palette$pch)
            }
        }
    } else if (type=="LOH") {
        if (is.null(res$input$vcf)) 
            stop("runAbsoluteCN was run without a VCF file.")
        if (is.null(ids)) ids <- 1:length(res$results)
        for (i in ids) {
            r <- res$results[[i]]$SNV.posterior$beta$posterior
            idx <-!r$ML.SOMATIC
            tmp <- cbind(levels(r[,1]), match(levels(r[,1]), r[idx, 1]))
            tmp <- tmp[complete.cases(tmp),,drop=FALSE]
            if (is.null(r)) next
            par(mfrow=c(3,1))
            main <- paste(
                "Purity:", round(res$results[[i]]$purity[[1]], digits=2), 
                " Tumor Ploidy:", round(res$results[[i]]$ploidy, digits=3),
                " SNV Log-Likelihood:", 
                    round(res$results[[i]]$SNV.posterior$beta$llik, digits=2),
                " Mean Coverage:", 
                    paste(round(apply(geno(res$input$vcf)$DP,2,mean)),
                    collapse=";") 
            )
            mycol <- ifelse(as.numeric(r[!r$ML.SOMATIC,1]) %% 2, "#E41A1C",
                "#377EB8")
            loh <- res$results[[i]]$SNV.posterior$beta$loh$output
            plot( do.call(c, lapply(1:nrow(loh), function(i) 
                rep(loh$seg.mean[i], loh$num.mark[i]))),
                ylab="Fraction LOH (segmented)", xlab="SNV Index",
                type="l", main=main, ...)
            axis(side=3, at=tmp[,2], labels=tmp[,1], tick=FALSE, padj=1)
            abline(v=tmp[,2], lty=3, col="grey")
            plot(r$ML.M[!r$ML.SOMATIC], 
                ylab="Maximum Likelihood SNV Multiplicity", 
                xlab="SNV Index", 
                ylim=c(0,min(7, max(r$ML.C[!r$ML.SOMATIC]))), col=mycol,...)
            abline(v=tmp[,2], lty=3, col="grey")
            if (sum(r$ML.SOMATIC)>0) {
                s <- r[!idx,paste("SOMATIC.M", 1:7, sep="")]
                colnames(s) <- gsub("SOMATIC.","",colnames(s))
                boxplot(s, xlab="Multiplicity somatic mutations", 
                    ylab="Posterior probability")
            }
        }
    } else if (type =="all") {
        plotAbs(res, type="overview")
        if (is.null(ids)) ids <- 1:length(res$results)
        for (i in ids) {
            par(mfrow=c(1,1))
            plotAbs(res,i, type="hist")
            if (!is.null(res$input$vcf)) {
                plotAbs(res,i, type="BAF",...)
                plotAbs(res,i, type="AF",...)
                plotAbs(res,i, type="LOH",...)
            }
        }    
    } else {
        mycol <-  ifelse(sapply(res$results, function(x) x$flag), "black", 
            "black")
        myfont <-  ifelse(sapply(res$results, function(x) x$flag), 3, 1)
        main <- NULL
    #    opar <- par(no.readonly=TRUE)
    #    par(mar= c(5, 4, 4, 4) + 0.1)
        #if ("darkgrey" %in% myfont) main="Italics: SCNA-fitting flagged."
        xc <- .matrixTotalPloidyToTumorPloidy(res$candidates$all)
        xc[is.infinite(xc)] <- min(xc[!is.infinite(xc)])
        xc[xc < quantile(xc, p=0.2)] <- quantile(xc, p=0.2)
        
        mycol.image <- colorRampPalette(rev(brewer.pal(n = 7, 
            name = "RdYlBu")))(100)
        image(as.numeric(colnames(xc)), as.numeric(rownames(xc)), 
            t(xc)-max(xc), col=mycol.image, xlab = "Purity", 
            ylab = "Ploidy",main = main,...)

        if (show.contour) contour(as.numeric(colnames(xc)), 
            as.numeric(rownames(xc)), t(xc), add=TRUE)

        text( sapply(res$results, function(x) min(0.9, x$purity))-0.02,
            sapply(res$results, function(x) x$ploidy)-0.1, 
            1:length(res$results), col=mycol, cex=2,font=myfont)
    #    par( opar ) 
    }
### Returns NULL
},ex=function() {
data(purecn.example.output)
plotAbs(purecn.example.output, type="overview")
# plot details for the maximum likelihood solution (rank 1)
plotAbs(purecn.example.output, 1, type="hist")
plotAbs(purecn.example.output, 1, type="BAF")
})    

.toLines <- function(
### "segments" already segmented log-ratios into a list for plotting
ss) {
    #get brakepoints
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

    cc <- do.call(cbind, lapply(1:length(tumor.ploidy.i), function(i) 
        ca[tumor.ploidy.i[[i]],i]))

    colnames(cc) <- colnames(ca)
    rownames(cc) <- rownames(ca)
    cc
}
