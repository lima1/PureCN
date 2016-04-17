createCurationFile <- structure(function(# Create file to curate PureCN results
### Function to create a CSV file that can be used to mark the correct solution
### in the output of a runAbsoluteCN() run.
file.rds,
### Output of the runAbsoluteCN() function, serialized with saveRDS()
overwrite.uncurated=TRUE
### Overwrite existing files unless flagged as "Curated".
) {
    res <- readRDS(file.rds)$results[[1]]

    d.f.curation <- data.frame(
        Sampleid=res$seg$ID[1], 
        Purity=res$purity, 
        Ploidy=res$ploidy, 
        Flagged=res$flag, 
        Failed=FALSE, 
        Curated=FALSE, 
        Comment=res$flag_comment
    )

    filename <- file.path(dirname(file.rds), 
        paste(gsub(".rds$", "", basename(file.rds)), "csv", sep="."))

    if (file.exists(filename)) {
        tmp <- read.csv(filename, as.is=TRUE)
        if (tmp$Curated[1]) {
            warning(filename, 
                " already exists and seems to be edited.",
                " Will not overwrite it.")
        } else if (!overwrite.uncurated) {
            warning(filename, " already exists. Will not overwrite it.")
        } else {
            write.csv(d.f.curation, file=filename, row.names=FALSE)
        }       
    } else {   
        write.csv(d.f.curation, file=filename, row.names=FALSE)
    }
    invisible(d.f.curation)
###A data.frame with the tumor purity and ploidy of the maximum likelihood
###solution 
},ex=function() {
data(purecn.example.output)
file.rds <- 'Sample1_PureCN.rds'
saveRDS(purecn.example.output, file=file.rds)
createCurationFile(file.rds) 
})
