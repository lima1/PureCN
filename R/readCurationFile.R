readCurationFile <- structure(function(#Read curation file
### Function that can be used to read the curated output of
### the runAbsoluteCN function.
file.rds,
### Output of the runAbsoluteCN() function, serialized with 
### saveRDS()
file.curation=gsub(".rds$", ".csv", file.rds),
### Filename of a curation file that points to the correct 
### tumor purity and ploidy solution.
remove.failed=FALSE,
### Do not return solutions that failed.
report.best.only=FALSE,
### Only return correct/best solution (useful on low memory
### machines when lots of samples are loaded)
min.ploidy=NULL,
### Minimum ploidy to be considered. If NULL, all. Can be 
### used to automatically ignore unlikely solutions.
max.ploidy=NULL
### Maximum ploidy to be considered. If NULL, all. Can be 
### used to automatically ignore unlikely solutions.
) {
    res <- readRDS(file.rds)
    curation <- read.csv(file.curation)
    if (curation$Failed) {
        if (remove.failed) return(NA)
        for (i in seq_along(res$results)) res$results[[i]]$failed <- TRUE
    } else {
        for (i in seq_along(res$results)) res$results[[i]]$failed <- FALSE
    }            
    diff.correct <- sapply(res$results, function(x) {
        abs(x$purity-curation$Purity) + (abs(x$ploidy-curation$Ploidy)/6)
    })
    idx.correct <- which.min(diff.correct)
    if (idx.correct != 1) {
        res$results[c(1,idx.correct)] <-  res$results[c(idx.correct, 1)]
    } 
    ploidy <- sapply(res$results, function(x) x$ploidy)
    if (is.null(min.ploidy)) min.ploidy <- min(ploidy)
    if (is.null(max.ploidy)) max.ploidy <- max(ploidy)

    idx <- which(ploidy>=min.ploidy & ploidy <= max.ploidy)
    res$results <- res$results[idx]
     
    if (report.best.only) {
        res$results <- res$results[1]
    }       
    res
### The return value of the corresponding runAbsoluteCN call, but with the
### results array manipulated according the curation CSV file and arguments
### of this function.
},ex=function() {
data(purecn.example.output)
file.rds <- 'Sample1_PureCN.rds'
createCurationFile(file.rds) 
# User can change the maximum likelihood solution manually in the generated 
# CSV file. The correct solution is then loaded with readCurationFile.
purecn.curated.example.output <-readCurationFile(file.rds) 
})    
