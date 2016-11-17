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
return.vcf=FALSE
### Returns an annotated \code{CollapsedVCF} object. Note that 
### this VCF will only contain variants not filtered out by the 
### \code{filterVcf} functions. Variants outside segments or intervals
### might be included or not depending on \code{\link{runAbsoluteCN}}
### arguments.
){
    pp   <- .addSymbols(res$results[[id]])
    if (return.vcf) {
        vcf <- res$input$vcf[
            res$results[[id]]$SNV.posterior$beta.model$vcf.ids]
        return(.annotatePosteriorsVcf(pp, vcf))
    }
    pp
### A \code{data.frame} or \code{CollapsedVCF} with SNV state 
### posterior probabilities.    
}, ex=function(){
data(purecn.example.output)
# the output data was created using a matched normal sample, but in case
# no matched normal is available, this will help predicting somatic vs. 
# germline status
purecn.snvs <- predictSomatic(purecn.example.output)

# write a VCF file
purecn.vcf <- predictSomatic(purecn.example.output, return.vcf=TRUE)
writeVcf(purecn.vcf, file="Sample1_PureCN.vcf")
})
    
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

    descriptionPp <- paste(ifelse(grepl("^S", colnames(pp)[idxColsPp]),
        "Somatic", "Germline"), gsub("^.*M", 
        " posterior probability, multiplicity ", 
            colnames(pp)[idxColsPp]), ".", sep="")

    descriptionPp <- gsub("GCONTHIGH", 
        " homozygous, reference allele contamination", descriptionPp)
    descriptionPp <- gsub("GCONTLOW", 
        " alt allele contamination", descriptionPp)
    descriptionPp <- gsub("GHOMOZYGOUS", 
        " homozygous", descriptionPp)

    newInfoPosterior <- DataFrame(
        Number=1, 
        Type="Float", 
        Description=descriptionPp, 
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
    idxColsMisc <- c("CN.SUBCLONAL", "CELLFRACTION",
    "Log.Ratio", "gene.symbol")
    newInfoMisc <- DataFrame(
        Number=c(0,1,1,1),
        Type=c("Flag", "Float", "String", "Float"),
        Description=c("Sub-clonal copy number gain", "Cellular fraction", 
            "Copy number log-ratio", "Gene Symbol" ),
        row.names=idxColsMisc
    )
    info(header(vcf)) <- rbind(info(header(vcf)), newInfoPosterior, newInfoMl,
newInfoMisc)
    idxColsMisc <- match(idxColsMisc, colnames(pp))
    pp[, idxColsPp] <- round( pp[, idxColsPp], digits=4)
    info(vcf) <- cbind(info(vcf), pp[, c(idxColsPp, idxCols, idxColsMisc)])
    vcf
}

