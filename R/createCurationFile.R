#' Create file to curate PureCN results
#' 
#' Function to create a CSV file that can be used to mark the correct solution
#' in the output of a \code{\link{runAbsoluteCN}} run.
#' 
#' 
#' @param file.rds Output of the \code{\link{runAbsoluteCN}} function,
#' serialized with \code{saveRDS}.
#' @param overwrite.uncurated Overwrite existing files unless flagged as
#' \sQuote{Curated}.
#' @param overwrite.curated Overwrite existing files even if flagged as
#' \sQuote{Curated}.
#' @return A \code{data.frame} with the tumor purity and ploidy of the maximum
#' likelihood solution.
#' @author Markus Riester
#' @seealso \code{\link{runAbsoluteCN}}
#' @examples
#' 
#' data(purecn.example.output)
#' file.rds <- "Sample1_PureCN.rds"
#' saveRDS(purecn.example.output, file=file.rds)
#' createCurationFile(file.rds) 
#' 
#' @export createCurationFile
#' @importFrom utils write.csv
createCurationFile <- function(file.rds, overwrite.uncurated = TRUE, 
    overwrite.curated=FALSE) {
    rds <- readRDS(file.rds)
    res <- rds$results[[1]]
    contamination <- res$SNV.posterior$posterior.contamination
    contamination <- if (is.null(contamination)) 0 else contamination
    d.f.curation <- data.frame(
        Sampleid=res$seg$ID[1], 
        Purity=res$purity, 
        Ploidy=res$ploidy, 
        Sex=.getSexFromRds(rds),
        Contamination=contamination,
        Flagged=res$flag, 
        Failed=FALSE, 
        Curated=FALSE, 
        Comment=res$flag_comment
    )

    filename <- file.path(dirname(file.rds), 
        paste(gsub(".rds$", "", basename(file.rds)), "csv", sep="."))

    if (file.exists(filename)) {
        tmp <- read.csv(filename, as.is=TRUE)
        if (tmp$Curated[1] && !overwrite.curated) {
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
}

.getSexFromRds <- function(rds) {
    if (!is.na(rds$input$sex) && !is.na(rds$input$sex.vcf)) {
        if (rds$input$sex == rds$input$sex.vcf) return(rds$input$sex)
        return(paste("Coverage:", rds$input$sex, "VCF:", rds$input$sex.vcf))     
    } 
    if (!is.na(rds$input$sex)) {
        return(rds$input$sex)
    }    
    return(rds$input$sex.vcf)
}
