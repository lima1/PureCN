#' Calling of amplifications and deletions
#' 
#' Function to extract major copy number alterations from a
#' \code{\link{runAbsoluteCN}} return object.
#' 
#' 
#' @param res Return object of the \code{\link{runAbsoluteCN}} function.
#' @param id Candidate solutions to be used. \code{id=1} will use the maximum
#' likelihood (or curated) solution.
#' @param cutoffs Copy numbers cutoffs to call losses, focal amplifications and
#' broad amplifications.
#' @param log.ratio.cutoffs Copy numbers log-ratio cutoffs to call losses and
#' amplifications in failed samples.
#' @param failed Indicates whether sample was failed. If \code{NULL}, use
#' available annotation, which can be set in the curation file.
#' @param all.genes If \code{FALSE}, then only return amplifications and
#' deletions passing the thresholds.
#' @return A \code{data.frame} with gene-level amplification and deletion
#' calls.
#' @author Markus Riester
#' @seealso \code{\link{runAbsoluteCN}}
#' @examples
#' 
#' data(purecn.example.output)
#' callAlterations(purecn.example.output)
#' callAlterations(purecn.example.output, all.genes=TRUE)["ESR2",]
#' 
#' @export callAlterations
callAlterations <- function(res, id = 1, cutoffs = c(0.5, 6, 7),
log.ratio.cutoffs = c(-0.9, 0.9), failed = NULL, all.genes = FALSE) {

    if (class(res$results[[id]]$gene.calls) != "data.frame") {
        .stopUserError("This function requires gene-level calls.\n",
            "Please add a column 'Gene' containing gene symbols to the ",
            "interval.file.")
    }
     
    amp.ids <- (res$results[[id]]$gene.calls$focal &
                res$results[[id]]$gene.calls$C >= cutoffs[2]) |
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

    bm <- res$results[[id]]$SNV.posterior
    if (!is.null(bm)) {
        segids <- bm$posteriors$seg.id
        calls$num.snps.segment <- sapply(calls$seg.id, function(i) 
            sum(segids==i,na.rm=TRUE))
        calls$M <- bm$posteriors$ML.M.SEGMENT[match(calls$seg.id, segids)] 
        calls$M.flagged <- bm$posteriors$M.SEGMENT.FLAGGED[match(calls$seg.id, segids)] 
        calls$loh <- bm$posteriors$ML.M.SEGMENT[match(calls$seg.id, segids)] == 0 
    }
    
    calls <- calls[, !grepl("^\\.",colnames(calls))]
       
    if (!all.genes) {
        return(calls[!is.na(calls$type),])
    }
    calls
}


#' Calling of amplifications and deletions from segmentations
#' 
#' This function can be used to obtain gene-level copy number calls from
#' segmentations. This is useful for comparing PureCN's segmentations with
#' segmentations obtained by different tools on the gene-level.  Segmentation
#' file can contain multiple samples.
#' 
#' 
#' @param sampleid The sampleid column in the segmentation file.
#' @param chr The chromosome column.
#' @param start The start positions of the segments.
#' @param end The end positions of the segments.
#' @param num.mark Optionally, the number of probes or markers in each segment.
#' @param seg.mean The segment mean.
#' @param C The segment integer copy number.
#' @param interval.file A mapping file that assigns GC content and gene symbols
#' to each exon in the coverage files. Used for generating gene-level calls.
#' First column in format CHR:START-END. Second column GC content (0 to 1).
#' Third column gene symbol. This file is generated with the 
#' \code{\link{preprocessIntervals}} function.
#' @param fun.focal Function for identifying focal amplifications. Defaults to
#' \code{\link{findFocal}}.
#' @param args.focal Arguments for focal amplification function.
#' @param \dots Arguments passed to \code{\link{callAlterations}}.
#' @return A list of \code{\link{callAlterations}} \code{data.frame} objects,
#' one for each sample.
#' @author Markus Riester
#' @examples
#' 
#' data(purecn.example.output)
#' seg <- purecn.example.output$results[[1]]$seg
#' interval.file <- system.file("extdata", "example_intervals.txt", 
#'         package = "PureCN")
#' 
#' calls <- callAlterationsFromSegmentation(sampleid=seg$ID, chr=seg$chrom,
#'     start=seg$loc.start, end=seg$loc.end, num.mark=seg$num.mark,
#'     seg.mean=seg$seg.mean, C=seg$C, interval.file=interval.file)
#' 
#' @export callAlterationsFromSegmentation
callAlterationsFromSegmentation <- function(sampleid, chr, start, end, 
    num.mark = NA, seg.mean, C, interval.file, fun.focal=findFocal,
    args.focal=list(), ...){
    seg <- data.frame(
        ID=sampleid,
        chrom=chr,
        loc.start=start,
        loc.end=end,
        num.mark=num.mark,
        seg.mean=seg.mean
    )    
    seg.adjusted <- data.frame(seg, C=C, size=seg$loc.end-seg$loc.start+1)

    tumor <- .addGCData(.gcGeneToCoverage(interval.file, 16), interval.file)
    chr.hash <- .getChrHash(seqlevels(tumor))
    segs <- split(seg.adjusted, seg$ID)
    gene.calls <- lapply(segs, function(s) {
        log.ratio <- .createFakeLogRatios(tumor, s[,1:6], s$ID[1], chr.hash)
        .getGeneCalls(s, tumor, log.ratio, fun.focal, args.focal, chr.hash)
    })
    res <- lapply(gene.calls, function(x) list(results=list(list(gene.calls=x, failed=FALSE))))
    lapply(res, callAlterations, ...)
}    
