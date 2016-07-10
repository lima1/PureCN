library('getopt')

spec <- matrix(c(
'help' , 'h', 0, "logical",
'normal', 'n', 1, "character",
'tumor', 't', 1, "character",
'vcf', 'v', 1, "character",
'genome' , 'g', 1, "character",
'gcgene' , 'c', 1, "character",
'segfile' , 'f', 1, "character",
'snpblacklist' , 's', 1, "character",
'exonweightfile' , 'e', 1, "character",
'normaldb' , 'd', 1, "character",
'outdir' , 'o', 1, "character",
'sampleid' , 'i', 1, "character"
), byrow=TRUE, ncol=4)
opt <- getopt(spec)

if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE))
    q(status=1)
}

gatk.normal.file <- opt$normal
gatk.tumor.file <- opt$tumor
tumor.vcf <- opt$vcf
genome <- opt$genome
gc.gene.file <- opt$gcgene
snp.blacklist <- opt$snpblacklist
seg.file <- opt$segfile
exon.weight.file <- opt$exonweightfile
normalDB <- opt$normaldb
sampleid <- opt$sampleid
outdir <- opt$outdir

PureCN <- function(
gatk.tumor.file, gatk.normal.file=NULL, tumor.vcf, genome,
gc.gene.file=NULL, seg.file=NULL, snp.blacklist=NULL, exon.weight.file=NULL, normalDB=NULL, 
sampleid, outdir
) {

    file.rds <- paste(outdir,"/",  sampleid, '_abs.rds', sep='')

    if (!is.null(normalDB)) {
        message("normalDB: ", normalDB)
        normalDB <- readRDS(normalDB)
        gatk.normal.file <- findBestNormal(gatk.tumor.file, normalDB)
    } else if (is.null(gatk.normal.file) && is.null(seg.file)) {
        stop("Need either normalDB or gatk.normal.file")
    }    
    message(paste('Best Normal:', gatk.normal.file))
    pdf(paste(outdir,"/", sampleid, '_abs_segmentation.pdf', sep=''), 
        width=10, height=12)

    ret <- runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
            gatk.tumor.file=gatk.tumor.file, vcf.file=tumor.vcf,
            sampleid=sampleid, gc.gene.file=gc.gene.file, plot.cnv=TRUE,
            genome=genome, seg.file=seg.file,
            args.filterVcf=list(snp.blacklist=snp.blacklist), 
            args.segmentation=list(exon.weight.file=exon.weight.file), 
            post.optimize=FALSE)
    dev.off()

    save(ret, file=paste(outdir,"/", sampleid, '_abs.rda', sep=''))
    saveRDS(ret, file=file.rds)
    createCurationFile(file.rds)
    pdf(paste(outdir,"/",sampleid, '_abs.pdf', sep=''), width=10, height=12)
    plotAbs(ret, type='all')
    dev.off()
}


outdir <- normalizePath(outdir, mustWork=TRUE)

if (is.null(seg.file)) {
    gatk.tumor.file <- normalizePath(gatk.tumor.file, mustWork=TRUE)
}

if (is.null(sampleid)) stop("Need sampleid.")

library(PureCN)

PureCN(gatk.tumor.file, gatk.normal.file, tumor.vcf, genome,
gc.gene.file, seg.file=seg.file, snp.blacklist=snp.blacklist, 
exon.weight.file=exon.weight.file, normalDB=normalDB, 
sampleid, outdir) 

