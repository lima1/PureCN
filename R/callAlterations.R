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
### Indicates whether sample was failed. If \code{NULL}, use available 
### annotation, which can be set in the curation file.
all.genes=FALSE,
### If \code{FALSE}, then only return amplifications and deletions 
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

    bm <- res$results[[id]]$SNV.posterior$beta.model
    if (!is.null(bm)) {
        segids <- bm$segment.ids
        calls$num.snps.segment <- sapply(calls$seg.id, function(i) 
            sum(segids==i,na.rm=TRUE))
        calls$loh <- bm$posterior$ML.M.Segment[match(calls$seg.id, segids)] == 0 
    }
    
    if (!all.genes) {
        return(calls[!is.na(calls$type),])
    }
    calls
### A \code{data.frame} with gene-level amplification and 
### deletion calls.
},ex=function() {
data(purecn.example.output)
callAlterations(purecn.example.output)
callAlterations(purecn.example.output, all.genes=TRUE)["ESR2",]
})        

callAlterationsFromSegmentation <- structure(
function(# Calling of amplifications and deletions from segmentations
### This function can be used to obtain gene-level copy number calls from
### segmentations. This is useful for comparing PureCN's segmentations with 
### segmentations obtained by different tools on the gene-level. 
### Segmentation file can contain multiple samples.
sampleid,
### The sampleid column in the segmentation file.
chr,
### The chromosome column.
start,
### The start positions of the segments.
end,
### The end positions of the segments.
num.mark=NA,
### Optionally, the number of probes or markers in each segment.
seg.mean,
### The segment mean.
C,
### The segment integer copy number.
gc.gene.file,
### A mapping file that assigns GC content and gene symbols 
### to each exon in the coverage files. Used for generating gene-level calls. 
### First column in format CHR:START-END. Second column GC content (0 to 1). 
### Third column gene symbol. This file
### can be generated with the \sQuote{GATK GCContentByInterval} tool or 
### with the \code{\link{calculateGCContentByInterval}} function.
fun.focal=findFocal,
### Function for identifying focal amplifications. Defaults to 
### \code{\link{findFocal}}.
args.focal=list(),
### Arguments for focal amplification function.
...
### Arguments passed to \code{\link{callAlterations}}.
){
    seg <- data.frame(
        ID=sampleid,
        chrom=chr,
        loc.start=start,
        loc.end=end,
        num.mark=num.mark,
        seg.mean=seg.mean
    )    
    seg.adjusted <- data.frame(seg, C=C, size=seg$loc.end-seg$loc.start+1)

    gc.data <- read.delim(gc.gene.file, as.is=TRUE)
    tumor <- .gcGeneToCoverage(gc.gene.file, 16)
    chr.hash <- .getChrHash(tumor$chr)
    segs <- split(seg.adjusted, seg$ID)
    gene.calls <- lapply(segs, function(s) {
        log.ratio <- .createFakeLogRatios(tumor, s[,1:6], chr.hash)
        .getGeneCalls(s, gc.data, log.ratio, fun.focal, args.focal, chr.hash)
    })
    res <- lapply(gene.calls, function(x) list(results=list(list(gene.calls=x, failed=FALSE))))
    lapply(res, callAlterations, ...)
### A list of \code{\link{callAlterations}} \code{data.frame} objects, 
### one for each sample.    
}, ex=function() {
    data(purecn.example.output)
    seg <- purecn.example.output$results[[1]]$seg
    gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
            package = "PureCN")

    calls <- callAlterationsFromSegmentation(sampleid=seg$ID, chr=seg$chrom,
        start=seg$loc.start, end=seg$loc.end, num.mark=seg$num.mark,
        seg.mean=seg$seg.mean, C=seg$C, gc.gene.file=gc.gene.file)
})
