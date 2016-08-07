segmentationPSCBS <-
structure(function(# PSCBS segmentation
### This function is deprecated and will be made defunct in the next
### Bioconductor release. Please use segmentationCBS instead.
normal, 
### GATK coverage data for normal sample.
tumor,  
### GATK coverage data for tumor sample.
log.ratio, 
### Copy number log-ratios, one for each exon in coverage file.
plot.cnv, 
### Segmentation plots.
min.coverage, 
### Minimum coverage in both normal and tumor.
sampleid=sampleid,
### Sample id, used in output files.
exon.weight.file=NULL,
### Can be used to assign weights to exons. NOT SUPPORTED YET.
flavor="tcn&dh",
### Flavor value for PSCBS. See segmentByNonPairedPSCBS.
tauA=0.03,
### tauA argument for PSCBS. See segmentByNonPairedPSCBS.
vcf=NULL,
### Optional VCF object with germline allelic ratios.
tumor.id.in.vcf=1,
### Id of tumor in case multiple samples are stored in VCF.
normal.id.in.vcf=NULL,
### Id of normal in in VCF. If NULL, use unpaired PSCBS.
max.segments=NULL,
### If not NULL, try a higher undo.SD parameter if number of
### segments exceeds the threshold.
chr.hash=NULL,
### Mapping of non-numerical chromsome names to numerical names
### (e.g. chr1 to 1, chr2 to 2, etc.). If NULL, assume chromsomes
### are properly ordered.   
verbose=TRUE,
### Verbose output.
...
### Additional parameters passed to the segmentByNonPairedPSCBS function.
) {
    debug <- TRUE
    .Deprecated("segmentationCBS")
        
    if (!requireNamespace("PSCBS", quietly = TRUE)) {
        .stopUserError("segmentationPSCBS requires the PSCBS package.")
    }

    if (is.null(chr.hash)) chr.hash <- .getChrHash(tumor$chr)

    exon.weights <- NULL
    if (!is.null(exon.weight.file)) {
        exon.weights <- read.delim(exon.weight.file, as.is=TRUE)
        exon.weights <- exon.weights[match(as.character(tumor[,1]), 
            exon.weights[,1]),2]
        if (verbose) message(
            "Exon weights found, but currently not supported by PSCBS.")
    }
    well.covered.exon.idx <- ((normal$average.coverage > min.coverage) &
        (tumor$average.coverage > min.coverage)) | 
        ((normal$average.coverage > 1.5 * min.coverage) &  
        (tumor$average.coverage > 0.5 * min.coverage))

    #MR: fix for missing chrX/Y 
    well.covered.exon.idx[is.na(well.covered.exon.idx)] <- FALSE
    tumor <- tumor[well.covered.exon.idx,]
    log.ratio <- log.ratio[well.covered.exon.idx]
    exon.gr <- GRanges(seqnames=tumor$chr, 
        IRanges(start=tumor$probe_start, end=tumor$probe_end))
    ov <- findOverlaps(vcf, exon.gr)
    d.f <- cbind(tumor[subjectHits(ov),], 
        CT=2^(log.ratio+1)[subjectHits(ov)], 
        betaT=unlist(geno(vcf[queryHits(ov)])$FA[,tumor.id.in.vcf]), 
        betaN=NA,
        x=start(vcf[queryHits(ov)]) )
    
    if (!is.null(normal.id.in.vcf)) {
        d.f$betaN <- unlist(geno(vcf[queryHits(ov)])$FA[,normal.id.in.vcf])
    }
         
    d.f.2 <- cbind(tumor[-subjectHits(ov),], 
        CT=2^(log.ratio+1)[-subjectHits(ov)], betaT=NA, betaN=NA,
        x=tumor$probe_start[-subjectHits(ov)] )
    
    d.f.3 <- rbind(d.f, d.f.2)
    d.f.3 <- d.f.3[order(.strip.chr.name(d.f.3$chr, chr.hash), d.f.3$x),]
    d.f <- d.f.3
    colnames(d.f)[2] <- "chromosome"
    d.f$chromosome <- .strip.chr.name(d.f$chromosome, chr.hash)
    #if (!is.null(normal.id.in.vcf)) {
    #    seg <- PSCBS::segmentByPairedPSCBS(d.f, tauA=tauA, 
    #        flavor=flavor, ...)
    #} else {
        seg <- PSCBS::segmentByNonPairedPSCBS(d.f, tauA=tauA, 
            flavor=flavor, ...)
    #}    
    if (plot.cnv) PSCBS::plotTracks(seg)
    .PSCBSoutput2DNAcopy(seg, sampleid)
### A list with elements seg and size. "seg" contains the 
### segmentation, "size" the size of all segments in base pairs.    
},ex=function() {
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")
vcf.file <- system.file("extdata", "example_vcf.vcf", 
    package="PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
    package="PureCN")

# The max.candidate.solutions, max.ploidy and test.purity parameters are set to
# non-default values to speed-up this example.  This is not a good idea for real
# samples.
ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
    gatk.tumor.file=gatk.tumor.file, vcf.file=vcf.file, genome="hg19",
    sampleid='Sample1', gc.gene.file=gc.gene.file,
    fun.segmentation=segmentationPSCBS, max.ploidy=4,
    test.purity=seq(0.3,0.7,by=0.05), max.candidate.solutions=1)
})    
    
.PSCBSoutput2DNAcopy <- function(seg, sampleid) {
    sx <- cbind(ID=sampleid, seg$output[!is.na(seg$output$tcnMean),])
    sx <- sx[,c("ID", "chromosome", "tcnStart", "tcnEnd", "tcnNbrOfLoci", 
        "tcnMean")]
    colnames(sx) <- c("ID", "chrom", "loc.start",  "loc.end", "num.mark", 
        "seg.mean")
    sx$seg.mean <- log2(sx$seg.mean/2)
    sx
}
