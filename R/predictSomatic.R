#' Predict germline vs. somatic status
#'
#' This function takes as input the output of a \code{\link{runAbsoluteCN}} run
#' and provides SNV posterior probabilities for all possible states.
#'
#'
#' @param res Return object of the \code{\link{runAbsoluteCN}} function.
#' @param id Candidate solutions to be analyzed. \code{id=1} will analyze the
#' maximum likelihood solution.
#' @param return.vcf Returns an annotated \code{CollapsedVCF} object. Note that
#' this VCF will only contain variants not filtered out by the \code{filterVcf}
#' functions. Variants outside segments or intervals might be included or not
#' depending on \code{\link{runAbsoluteCN}} arguments.
#' @return A \code{data.frame} or \code{CollapsedVCF} with SNV state posterior
#' probabilities.
#' @author Markus Riester
#' @seealso \code{\link{runAbsoluteCN}}
#' @examples
#'
#' data(purecn.example.output)
#' # the output data was created using a matched normal sample, but in case
#' # no matched normal is available, this will help predicting somatic vs. 
#' # germline status
#' purecnSnvs <- predictSomatic(purecn.example.output)
#'
#' # Prefer GRanges?
#' purecnSnvs <- GRanges(predictSomatic(purecn.example.output))
#' 
#' # write a VCF file
#' purecnVcf <- predictSomatic(purecn.example.output, return.vcf=TRUE)
#' writeVcf(purecnVcf, file = "Sample1_PureCN.vcf")
#' 
#' @export predictSomatic
predictSomatic <- function(res, id = 1, return.vcf = FALSE) {
    pp   <- .addSymbols(res$results[[id]])
    if (return.vcf) {
        vcf <- res$input$vcf[
            res$results[[id]]$SNV.posterior$vcf.ids]
        return(.annotatePosteriorsVcf(pp, vcf))
    }
    pp
}

.calcMultisamplePosteriors <- function(ret.list) {
    pool <- Reduce(union, lapply(ret.list, function(r)
        as.character(rowRanges(r$input$vcf))))
    pos <- lapply(ret.list, function(r) match(as.character(GRanges(r$results[[1]]$SNV.posterior$posteriors)), pool))
    pos.align <- lapply(seq_along(pool), function(i)
                           lapply(pos, function(j) which(j == i)))
    likelihoods <- lapply(pos.align, function(i) do.call(rbind, lapply(seq_along(i), function(j)
              ret.list[[j]]$results[[1]]$SNV.posterior$likelihoods[i[[j]],])))
    idx <- grep("SOMATIC.M", colnames(likelihoods[[1]]))
    pp <- sapply(likelihoods, function(l) 
                 sum(apply(l,1,function(x) sum(x[idx]))) /
                 sum(apply(l,1,function(x) sum(x))))
    names(pp) <- pool
    lapply(ret.list, function(x) {
        ps <- predictSomatic(x)
        keys <- as.character(GRanges(ps))
        data.frame(Sampleid = x$input$sampleid,
                   ps,
                   MULTI.POSTERIOR.SOMATIC = pp[keys])
    })
}       

.addSymbols <- function(result) {
    if (is(result$gene.calls, "data.frame")) {
        g.gr <- GRanges(result$gene.calls)
        p.gr <- GRanges(result$SNV.posterior$posteriors)
        ov <- findOverlaps(p.gr, g.gr)
        
        result$SNV.posterior$posteriors$gene.symbol <- NA
        result$SNV.posterior$posteriors$gene.symbol[queryHits(ov)] <- 
            rownames(result$gene.calls)[subjectHits(ov)]
    }    
    result$SNV.posterior$posteriors
}
    
.annotatePosteriorsVcf <- function(pp, vcf) {
    if (nrow(pp) != nrow(vcf)) {
         .stopRuntimeError("Posteriors and filtered VCF do not align.")
    }
    if (is.null(pp$gene.symbol)) pp$gene.symbol <- "."

    idxColsPp <- grep("^SOMATIC|^GERMLINE", colnames(pp))
    colnames(pp)[idxColsPp] <- gsub("OMATIC\\.|ERMLINE\\.", "",
        colnames(pp)[idxColsPp])
    colnames(pp)[colnames(pp) == "CN.SUBCLONAL"] <- "CS"
    colnames(pp)[colnames(pp) == "CELLFRACTION"] <- "CF"
    colnames(pp)[colnames(pp) == "POSTERIOR.SOMATIC"] <- "PS"
    colnames(pp)[colnames(pp) == "log.ratio"] <- "LR"
    colnames(pp)[colnames(pp) == "gene.symbol"] <- "GS"

    descriptionPp <- paste0(ifelse(grepl("^S", colnames(pp)[idxColsPp]),
        "Somatic", "Germline"), gsub("^.*M(\\d)", 
        " posterior probability, multiplicity \\1", 
            colnames(pp)[idxColsPp]))

    descriptionPp <- gsub("GCONTHIGH", 
        " homozygous, reference allele contamination", descriptionPp)
    descriptionPp <- gsub("GCONTLOW", 
        " alt allele contamination", descriptionPp)
    descriptionPp <- gsub("GHOMOZYGOUS", 
        " homozygous", descriptionPp)
    idxCols <- grep("^ML", colnames(pp))
    prefix <- .getPureCNPrefixVcf(vcf)
    colnames(pp) <- paste0(prefix, colnames(pp))
    newInfoPosterior <- DataFrame(
        Number = 1, 
        Type = "Float", 
        Description = descriptionPp, 
        row.names = colnames(pp)[idxColsPp]
    )
    infoType <- as.character(sapply(sapply(pp[idxCols], class), function(x)
ifelse(x=="logical", "Flag", ifelse(x=="integer", "Integer", "Float"))))
    newInfoMl <- DataFrame(
        Number = ifelse(infoType == "Flag", 0, 1),
        Type = infoType, 
        Description = c(
        "Maximum likelihood state is a somatic state",
        "Maximum likelihood multiplicity",
        "Maximum likelihood integer copy number",
        "Maximum likelihood minor segment copy number",
        "Expected allelic fraction of the maximum likelihood state",
        "TRUE if segment is most likely in LOH"
        ),
        row.names = colnames(pp)[idxCols]
    )

    idxColsMisc <- paste0(prefix, c("CS", "CF", "PS", "LR", "GS", "FLAGGED"))
    newInfoMisc <- DataFrame(
        Number = c(0, 1, 1, 1, 1, 0),
        Type = c("Flag", "Float", "Float", "Float", "String", "Flag"),
        Description = c("Sub-clonal copy number gain", "Cellular fraction",
            "Posterior probability variant is somatic mutation.",
            "Copy number log-ratio", "Gene Symbol", "QC Flagged"),
        row.names = idxColsMisc
    )
    info(header(vcf)) <- rbind(info(header(vcf)), newInfoPosterior, newInfoMl,
        newInfoMisc)
    idxColsMisc <- match(idxColsMisc, colnames(pp))
    pp[, idxColsPp] <- round(pp[, idxColsPp], digits = 4)
    pp[, idxColsMisc[2:4]] <- round(pp[, idxColsMisc[2:4]], digits = 4)
  
    info(vcf) <- cbind(info(vcf), pp[, c(idxColsPp, idxCols, idxColsMisc)])
    vcf
}
