predictSomatic <-
structure(function(#Predict germline vs. somatic status
### This function takes as input the output of a 
### runAbsoluteCN run and annotates variants more correctly as germline 
### vs. somatic by inferring maternal and paternal chromosome numbers.
res, 
### Return object of the runAbsoluteCN() function.
id=1, 
### Candidate solutions to be analyzed. id=1 will analyze the 
### maximum likelihood solution.
cutoff=0.1
### Exclude maternal/paternal chromosome number in segment if 
### posterior probability is lower than cutoff.
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
    pp.new
### A data.frame with adjusted SNV state posterior probabilities.    
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
    

