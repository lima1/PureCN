callAlterations <- structure(
function(# Calling of amplifications and deletions
### Function to extract major copy number alterations from a 
### \code{\link{runAbsoluteCN}} return object.
res,
### Return object of the \code{\link{runAbsoluteCN}} function.
##seealso<< \code{\link{runAbsoluteCN}}
id=1,
### Candidate solutions to be used. \code{id=1} will use the 
### maximum likelihood (or curated) solution.
cutoffs=c(0.5,6,7),
### Copy numbers cutoffs to call losses, focal amplifications 
### and broad amplifications.
log.ratio.cutoffs=c(-0.9,0.9),
### Copy numbers log-ratio cutoffs to call losses and amplifications 
### in failed samples.
failed=NULL,
### Indicates whether sample was failed. If NULL, use available 
### annotation, which can be set in the curation file.
all.genes=FALSE,
### If FALSE, then only return amplifications and deletions 
### passing the thresholds.
...) {

    if (class(res$results[[id]]$gene.calls) != "data.frame") {
        .stopUserError("This function requires gene-level calls.\n",
            "Please add a column 'Gene' containing gene symbols to the ",
            "gc.gene.file.")
    }
     
    amp.ids <- ( res$results[[id]]$gene.calls$focal & 
                 res$results[[id]]$gene.calls$C >= cutoffs[2] ) |
                 res$results[[id]]$gene.calls$C >= cutoffs[3] 

    del.ids <- res$results[[id]]$gene.calls$C < cutoffs[1]

    if (is.null(failed)) failed <- res$results[[id]]$failed

    if (failed) {
        amp.ids <- res$results[[id]]$gene.calls$gene.mean >= 
            log.ratio.cutoffs[2] 
        del.ids <- res$results[[id]]$gene.calls$gene.mean < 
            log.ratio.cutoffs[1]
    }

    calls <- res$results[[1]]$gene.calls
    calls$type <- NA
    calls$type[amp.ids] <- "AMPLIFICATION"
    calls$type[del.ids] <- "DELETION"
    
    if (!all.genes) {
        return(calls[!is.na(calls$type),])
    }
    calls
### A data.frame with gene-level amplification and deletion calls.
},ex=function() {
data(purecn.example.output)
callAlterations(purecn.example.output)
callAlterations(purecn.example.output, all.genes=TRUE)["ESR2",]
})        
