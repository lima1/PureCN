calculateGCContentByInterval <- structure(
function(# Calculates GC content by interval
### Uses \code{scanFa} from the Rsamtools package to retrieve GC 
### content of intervals in a reference FASTA file.
interval.file,
### File specifying the intervals. Interval is expected in 
### first column in format CHR:START-END. BED files with 
### coordinates in the first three columns are supported as well.
reference.file,
### Reference FASTA file.
output.file=NULL,
### Optionally, write GC content file. 
...
### Additional parameters passed to the \code{read.delim}
### function that reads the \code{interval.file}.
) {
    interval <- read.delim(interval.file, as.is=TRUE, ...)
    isBedFile <- .checkBedFile(interval)
    if (isBedFile) {
       interval$Target <- paste0(interval[,1],":", interval[,2],"-", 
        interval[,3])
    } else {
        colnames(interval)[1] <- "Target"
    }
    pos <- as.data.frame(do.call(rbind, strsplit(interval$Target, ":|-")), 
        stringsAsFactors = FALSE)
    interval.gr <- GRanges(seqnames = pos[,1], 
        IRanges(start = as.numeric(pos[,2]), end = as.numeric(pos[,3])))
    if (min(start(interval.gr)) < 1) {
        .stopUserError("Interval coordinates should start at 1, not at 0.")
    }    
    x <- scanFa(reference.file, interval.gr)
    GC.count <- letterFrequency(x,"GC")
    all.count <- letterFrequency(x,"ATGC")
    gc <- data.frame(
        Target=interval$Target,
        gc_bias=as.vector(ifelse(all.count==0,NA,GC.count/all.count))
    )
    if (!is.null(output.file)) {
        write.table(gc, file=output.file, row.names=FALSE, quote=FALSE, 
            sep="\t")
    }    
    invisible(gc)
### Returns GC content by interval.
}, ex=function() {
reference.file <- system.file("extdata", "ex2_reference.fa", 
    package="PureCN", mustWork = TRUE)
interval.file <- system.file("extdata", "ex2_intervals.txt", 
    package="PureCN", mustWork = TRUE)
calculateGCContentByInterval(interval.file, reference.file, 
    output.file="gc_file.txt")
}) 


.checkBedFile <- function(interval) {
    seqStyle <- try(seqlevelsStyle(interval[,1]), silent=TRUE)
    return(class(seqStyle) != "try-error") 
}
