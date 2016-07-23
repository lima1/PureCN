predictSomatic <-
structure(function(#Predict germline vs. somatic status
### This function takes as input the output of a 
### \code{\link{runAbsoluteCN}} run and provides SNV posterior probabilities
### for all possible states. 
res, 
### Return object of the \code{\link{runAbsoluteCN}} function.
##seealso<< \code{\link{runAbsoluteCN}}
id=1, 
### Candidate solutions to be analyzed. \code{id=1} will analyze the 
### maximum likelihood solution.
cutoff=0.1,
### Exclude maternal/paternal chromosome number in segment if 
### posterior probability is lower than cutoff.
return.vcf=FALSE
### Returns an annotated \code{CollapsedVCF} object. Note that 
### this VCF will only contain variants not filtered out by the 
### \code{filterVcf} functions. Variants outside segments or intervals
### might be included or not depending on \code{\link{runAbsoluteCN}}
### arguments.
){
    llik <- res$results[[id]]$SNV.posterior$beta.model$likelihoods
    pp   <- .addSymbols(res$results[[id]])
    purity <- res$results[[id]]$purity

    llik.segment.sums <- .calcMpriorGermline(
        res$results[[id]]$SNV.posterior$beta.model)

    llik.segment.delete <- lapply(llik.segment.sums, function(x) x < cutoff)
    germline.cols.pp <- grep("GERMLINE.M",colnames(pp))
    pp.new <- pp
    for (i in names(llik.segment.delete)) {
         row.ids <- which(
            res$results[[id]]$SNV.posterior$beta.model$segment.ids == i)
         pp.new[row.ids, germline.cols.pp[llik.segment.delete[[i]]]] <- 0
    }
    min.ar <- 0
    if (sum(pp.new$ML.SOMATIC)>10) {
        min.ar <- min(pp.new$ML.AR[pp.new$ML.SOMATIC])
    }
    # set everything that does not match well and has low allelic fraction
    # as somatic 
    pp.new[pp$AR < min.ar & rowSums(llik)<1,"SOMATIC.M1"] <- 1

    snv.posteriors <-  pp.new[,4:21] 
    snv.posteriors <- snv.posteriors/rowSums(snv.posteriors) 
    xx <- .extractMLSNVState(snv.posteriors)
    pp.new[,4:21] <- snv.posteriors 
    pp.new[, colnames(xx)] <- xx
    pp.new$Cellfraction <-  (pp.new$AR/pp.new$ML.M)*
        (purity*pp.new$ML.C+2*(1-purity))/purity
    pp.new$Cellfraction[!pp.new$ML.SOMATIC] <- NA
    if (return.vcf) {
        vcf <- res$input$vcf[
            res$results[[id]]$SNV.posterior$beta.model$vcf.ids]
        return(.annotatePosteriorsVcf(pp.new, vcf))
    }
    pp.new
### A \code{data.frame} or \code{CollapsedVCF} with SNV state 
### posterior probabilities.    
}, ex=function(){
data(purecn.example.output)
# the output data was created using a matched normal sample, but in case
# no matched normal is available, this will help predicting somatic vs. 
# germline status
purecn.snvs <- predictSomatic(purecn.example.output)
})
    
.calcMpriorGermline <- function(model) {
    llik <- model$likelihoods
    llik.segment <- split(as.data.frame(llik), model$segment.ids)
    germline.cols <- grep("GERMLINE.M",colnames(llik))
    lapply(llik.segment, function(x) colSums(x[, germline.cols])/
        max(colSums(x[, germline.cols])))
}

.addSymbols <- function(result) {
    if (class(result$gene.calls) == "data.frame") {
        
    g.gr <- GRanges(seqnames=result$gene.calls$chr,
              IRanges(start=result$gene.calls$start,
                      end=result$gene.calls$end))
       
    p.gr <- GRanges(seqnames=result$SNV.posterior$beta.model$posteriors$chr,
              IRanges(start=result$SNV.posterior$beta.model$posteriors$start,
                      end=result$SNV.posterior$beta.model$posteriors$end))
    
    ov <- findOverlaps(p.gr, g.gr)
    
    result$SNV.posterior$beta.model$posteriors$gene.symbol <- NA
    result$SNV.posterior$beta.model$posteriors$gene.symbol[queryHits(ov)] <- 
        rownames(result$gene.calls)[subjectHits(ov)]
    }    
    result$SNV.posterior$beta.model$posteriors
}
    
.annotatePosteriorsVcf <- function(pp, vcf) {
    if (nrow(pp) != nrow(vcf)) {
         .stopRuntimeError("Posteriors and filtered VCF do not align.")
    }
    idxColsPp <- grep("^SOMATIC|^GERMLINE", colnames(pp))
    colnames(pp)[idxColsPp] <- gsub("OMATIC\\.|ERMLINE\\.", "",
        colnames(pp)[idxColsPp])

    newInfoPosterior <- DataFrame(
        Number=1, 
        Type="Float", 
        Description=gsub("^.*M", " posterior probability, multiplicity ", 
            colnames(pp)[idxColsPp]),
        row.names=colnames(pp)[idxColsPp]
    )
    idxCols <- grep("^ML", colnames(pp))
    infoType <- as.character(sapply(sapply(pp[idxCols], class), function(x)
ifelse(x=="logical", "Flag", ifelse(x=="integer", "Integer", "Float"))))
    newInfoMl <- DataFrame(
        Number=ifelse(infoType=="Flag", 0, 1),
        Type=infoType, 
        Description=gsub("^ML\\.", "Maximum likelihood estimate ", 
            colnames(pp)[idxCols]),
        row.names=colnames(pp)[idxCols]
    )
    idxColsMisc <- c("CN.Subclonal", "Log.Ratio", "gene.symbol", 
        "Cellfraction")
    newInfoMisc <- DataFrame(
        Number=c(0,1,1,1),
        Type=c("Flag", "Float", "String", "Float"),
        Description=c("Sub-clonal copy number gain", "Copy number log-ratio", 
            "Gene Symbol", "Cellular fraction"),
        row.names=idxColsMisc
    )
    info(header(vcf)) <- rbind(info(header(vcf)), newInfoPosterior, newInfoMl,
newInfoMisc)
    idxColsMisc <- match(idxColsMisc, colnames(pp))
    pp[, idxColsPp] <- round( pp[, idxColsPp], digits=4)
    info(vcf) <- cbind(info(vcf), pp[, c(idxColsPp, idxCols, idxColsMisc)])
    vcf
}

