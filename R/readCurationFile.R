#' Read curation file
#' 
#' Function that can be used to read the curated output of the
#' \code{\link{runAbsoluteCN}} function.
#' 
#' 
#' @param file.rds Output of the \code{\link{runAbsoluteCN}} function,
#' serialized with \code{saveRDS}.
#' @param file.curation Filename of a curation file that points to the correct
#' tumor purity and ploidy solution.
#' @param remove.failed Do not return solutions that failed.
#' @param report.best.only Only return correct/best solution (useful on low
#' memory machines when lots of samples are loaded).
#' @param min.ploidy Minimum ploidy to be considered. If \code{NULL}, all. Can
#' be used to automatically ignore unlikely solutions.
#' @param max.ploidy Maximum ploidy to be considered. If \code{NULL}, all. Can
#' be used to automatically ignore unlikely solutions.
#' @return The return value of the corresponding \code{\link{runAbsoluteCN}}
#' call, but with the results array manipulated according the curation CSV file
#' and arguments of this function.
#' @author Markus Riester
#' @seealso \code{\link{runAbsoluteCN} \link{createCurationFile}}
#' @examples
#' 
#' data(purecn.example.output)
#' file.rds <- "Sample1_PureCN.rds"
#' createCurationFile(file.rds) 
#' # User can change the maximum likelihood solution manually in the generated 
#' # CSV file. The correct solution is then loaded with readCurationFile.
#' purecn.curated.example.output <-readCurationFile(file.rds) 
#' 
#' @export readCurationFile
#' @importFrom utils read.csv
readCurationFile <- function(file.rds,
file.curation = gsub(".rds$", ".csv", file.rds),
remove.failed = FALSE, report.best.only=FALSE, min.ploidy = NULL,
max.ploidy = NULL) {
    flog.info("Reading %s...", file.rds)
    res <- readRDS(file.rds)
    if (!file.exists(file.curation)) {
        flog.warn("Curation file %s does not exist, creating one.", file.curation)
        createCurationFile(file.rds)
    }
    curation <- read.csv(file.curation, as.is=TRUE, nrows=1)
    .checkLogical <- function(field) {
        if (!is.logical(curation[[field]])) {
            .stopUserError("'", field, "' column in ", file.curation, 
                " not logical(1).")
        }
    }
    .checkLogical("Failed")
    .checkLogical("Curated")
    .checkLogical("Flagged")

    ## Mark all solutions as failed if sample is curated as failed
    if (curation$Failed) {
        if (remove.failed) return(NA)
        for (i in seq_along(res$results)) res$results[[i]]$failed <- TRUE
    } else {
        for (i in seq_along(res$results)) res$results[[i]]$failed <- FALSE
    }
    
    # Make sure purity and ploidy are numeric. Stop if not, not warn.
    curation$Purity <- suppressWarnings(as.numeric(curation$Purity))
    curation$Ploidy <- suppressWarnings(as.numeric(curation$Ploidy))
    
    if (is.na(curation$Purity) || is.na(curation$Ploidy) ||
        curation$Purity < 0 || curation$Purity > 1 ||
        curation$Ploidy < 0 || curation$Ploidy > 8) {
        .stopUserError("Purity or Ploidy not numeric or in expected range.")
    }    
    # Find purity/ploidy solution most similar to curation
    diffCurated <- vapply(res$results, function(x) {
        abs(x$purity-curation$Purity) + (abs(x$ploidy-curation$Ploidy)/6)
    }, double(1))
    idxCurated <- which.min(diffCurated)
    if (idxCurated != 1) {
        res$results[c(1,idxCurated)] <-  res$results[c(idxCurated, 1)]
    }
    
    ## Filter by ploidy if necessary
    ploidy <- sapply(res$results, function(x) x$ploidy)
    if (is.null(min.ploidy)) min.ploidy <- min(ploidy)
    if (is.null(max.ploidy)) max.ploidy <- max(ploidy)
    idxPloidyOk <- which(ploidy>=min.ploidy & ploidy <= max.ploidy)
    res$results <- res$results[idxPloidyOk]
     
    if (report.best.only) {
        res$results <- res$results[1]
    }
    res
}    
