#' Get regions of LOH
#' 
#' This function provides detailed LOH information by region.
#' 
#' 
#' @param res Return object of the \code{\link{runAbsoluteCN}} function.
#' @param id Candidate solution to extract LOH from. \code{id=1} will use the
#' maximum likelihood solution.
#' @param arm.cutoff Min fraction LOH on a chromosome arm to call whole arm
#' events.
#' @return Returns \code{data.frame} with LOH regions.
#' @author Markus Riester
#' @seealso \code{\link{runAbsoluteCN}}
#' @examples
#' 
#' data(purecn.example.output)
#' head(callLOH(purecn.example.output))
#' 
#' @export callLOH
callLOH <- function(res, id = 1, arm.cutoff = 0.9) {
    if (is.null(res$input$vcf)) {
        .stopUserError("runAbsoluteCN was run without a VCF file.")
    }
    chr.hash <- res$input$chr.hash
    armLocations <- .getArmLocations(res)
    if (!nrow(armLocations)) {
        .stopUserError("Centromere positions not available or matching.")
    }
    armLocationsGR <- GRanges(armLocations)
    seg <- res$results[[id]]$seg

    minorChrNumber <- res$results[[id]]$SNV.posterior$posteriors[, 
        c("seg.id", "ML.M.SEGMENT")]
    minorChrNumber <- minorChrNumber[!duplicated(minorChrNumber[,1]),]
    seg$M <- NA
    seg$M[minorChrNumber[,1]] <- minorChrNumber[,2]
    seg <- seg[!is.na(seg$M),]
    seg$chrom <- .add.chr.name(seg$chrom, chr.hash)
    
    # merge consecutive segments if they have same genotype
    i <- 1
    while (i < nrow(seg)) {
        key <- paste(seg$chrom, seg$C, seg$M)
        if (key[i] == key[i+1]) {
            seg$loc.end[i] <- seg$loc.end[i+1]
            seg$num.mark[i] <- seg$num.mark[i] + seg$num.mark[i+1]
            seg <- seg[-(i+1),]
            next
        }
        i <- i + 1
    }

    segGR <- GRanges(seg)

    ov <- findOverlaps(armLocationsGR, segGR)
    segLOH <-  cbind(seg[subjectHits(ov),], armLocations[queryHits(ov),])
    segLOH$loc.start <- apply(segLOH[,c("loc.start", "start")],1,max)
    segLOH$loc.end <- apply(segLOH[,c("loc.end", "end")],1,min)
    segLOH$size <- segLOH$loc.end-segLOH$loc.start+1
    segLOH$fraction.arm <- round(segLOH$size/armLocations$size[
        match(paste(segLOH$chrom, segLOH$arm),
        paste(armLocations$chrom, armLocations$arm))], digits=2)
    segLOH$type <- ""
    segLOH$type[segLOH$C == 2 & segLOH$M == 0] <- "COPY-NEUTRAL LOH"
    segLOH$type[segLOH$type=="" & segLOH$M == 0] <- "LOH"
    idx <- segLOH$fraction.arm > arm.cutoff & segLOH$type != ""
    segLOH$type[idx] <- paste("WHOLE ARM",
        segLOH$type)[idx]
    rownames(segLOH) <- NULL
    segLOH <- segLOH[, c("chrom", "loc.start", "loc.end", "arm", "C", "M", "type")]
    # standardize colnames
    colnames(segLOH)[1:3] <- c("chr", "start", "end")
    segLOH
}

.getCentromeres <- function(res) {
    # TODO remove this support for old data.frame centromeres in PureCN 1.12
    if (class(res$input$centromeres) == "GRanges" || 
            is.null(res$input$centromeres)) {
        return(res$input$centromeres)
    }    
    GRanges(res$input$centromeres)    
}
    
.getArmLocations <- function(res) {
    chr.hash <- res$input$chr.hash
    centromeres <- .getCentromeres(res)

    chromCoords <- suppressWarnings(t(vapply(split(
        start(res$input$vcf),
        as.character(seqnames(res$input$vcf))), function(x)
        c(min(x), max(x)), c(min=double(1), max=double(1)))))
    if (!is.null(centromeres)) {
        # split segments by centromere if available
        chromCoords <- chromCoords[as.integer(match(seqnames(centromeres), rownames(chromCoords))),]
        rownames(chromCoords) <- NULL

        centromeres <- cbind(data.frame(centromeres), chromCoords)
        centromeres <- centromeres[complete.cases(centromeres),]

        pArms <- centromeres[centromeres$min < centromeres$start,
            c("seqnames", "min", "start")]
        qArms <- centromeres[centromeres$max > centromeres$end,
            c("seqnames", "end", "max")]

        colnames(pArms) <- c("chrom", "start", "end")
        colnames(qArms) <- colnames(pArms)
        pArms$arm <- "p"
        if (nrow(qArms) == 0) {
          armLocations = pArms
        } else {
          qArms$arm <- "q"
          armLocations <- rbind(pArms, qArms)
        }

    } else {
        armLocations <- data.frame(
                            chrom=rownames(chromCoords),
                            start=chromCoords[,1],
                            end=chromCoords[,2],
                            arm="")
    }    
    armLocations <- armLocations[order(match(armLocations$chrom,
        chr.hash[,1])),]
    armLocations$size <- armLocations$end-armLocations$start+1
    armLocations
}
