library('getopt')

spec <- matrix(c(
'help',           'h', 0, "logical",
'force' ,         'f', 0, "logical",
'normal',         'n', 1, "character",
'tumor',          't', 1, "character",
'vcf',            'v', 1, "character",
'genome',         'g', 1, "character",
'gcgene',         'c', 1, "character",
'segfile',        'l', 1, "character",
'snpblacklist',   's', 1, "character",
'normal_panel',   'p', 1, "character",
'statsfile',      'a', 1, "character",
'targetweightfile', 'e', 1, "character",
'normaldb',       'd', 1, "character",
'pool',           'm', 1, "integer",
'postoptimize',   'z', 0, "logical",
'modelhomozygous','y', 0, "logical",
'outdir',         'o', 1, "character",
'outvcf',         'u', 1, "character",
'sampleid',       'i', 1, "character"
), byrow=TRUE, ncol=4)
opt <- getopt(spec)

if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE))
    q(status=1)
}

force <- !is.null(opt$force)
post.optimize <- !is.null(opt$postoptimize)
normal.coverage.file <- opt$normal
tumor.coverage.file <- opt$tumor
tumor.vcf <- opt$vcf
genome <- opt$genome
gc.gene.file <- opt$gcgene
snp.blacklist <- opt$snpblacklist
stats.file <- opt$statsfile
seg.file <- opt$segfile
target.weight.file <- opt$targetweightfile
normal.panel.vcf.file <- opt$normal_panel
normalDB <- opt$normaldb
sampleid <- opt$sampleid
outdir <- opt$outdir
outvcf <- opt$outvcf
pool <- opt$pool
model.homozygous <- opt$modelhomozygous

PureCN <- function(
tumor.coverage.file, normal.coverage.file=NULL, tumor.vcf, genome,
gc.gene.file=NULL, seg.file=NULL, snp.blacklist=NULL, stats.file=NULL,
target.weight.file=NULL, normal.panel.vcf.file=NULL, normalDB=NULL, sampleid,
post.optimize=FALSE, pool, model.homozygous, outvcf, outdir) {

    file.rds <- file.path(outdir, paste(sampleid, '_abs.rds', sep=''))

    if (file.exists(file.rds) && !force) {
        message(file.rds, 
            " already exists. Skipping... (--force will overwrite)")
        ret <- readRDS(file.rds)
    } else {    

        if (!is.null(normalDB)) {
            message("normalDB: ", normalDB)
            normalDB <- readRDS(normalDB)
            if (is.null(normal.coverage.file)) {
                if (!is.null(pool)) {
                    num.normals <- pool
                    pool <- TRUE
                } else {
                    num.normals <- 1
                    pool <- FALSE
                }    
                normal.coverage.file <- findBestNormal(tumor.coverage.file, 
                    normalDB, pool=pool, num.normals=num.normals)
                if (!pool) message(paste('Best Normal:', normal.coverage.file))
            }
        } else if (is.null(normal.coverage.file) && is.null(seg.file)) {
            stop("Need either normalDB or normal.coverage.file")
        }    

        pdf(paste(outdir,"/", sampleid, '_abs_segmentation.pdf', sep=''), 
            width=10, height=11)
        ret <- runAbsoluteCN(normal.coverage.file=normal.coverage.file, 
                tumor.coverage.file=tumor.coverage.file, vcf.file=tumor.vcf,
                sampleid=sampleid, gc.gene.file=gc.gene.file, plot.cnv=TRUE,
                genome=genome, seg.file=seg.file,
                args.filterVcf=list(snp.blacklist=snp.blacklist, 
                    stats.file=stats.file), 
                args.segmentation=list(target.weight.file=target.weight.file), 
                args.setMappingBiasVcf=
                    list(normal.panel.vcf.file=normal.panel.vcf.file),
                args.filterTargets=list(normalDB=normalDB),
                model.homozygous=model.homozygous,
                post.optimize=post.optimize)
        dev.off()
        saveRDS(ret, file=file.rds)
    }
    createCurationFile(file.rds)
    file.pdf <- file.path(outdir, paste(sampleid, '_abs.pdf', sep=''))
    pdf(file.pdf, width=10, height=11)
    plotAbs(ret, type='all')
    dev.off()
    if (!is.null(outvcf)) {
        if (basename(outvcf) != outvcf) {
            warning(outvcf, " contains path. Will not put it in ", outdir)
        } else {
            outvcf <- file.path(outdir, outvcf)
        }    
        vcfanno <- predictSomatic(ret, return.vcf=TRUE, 
            vcf.field.prefix="PureCN_")
        writeVcf(vcfanno, file=outvcf)    
    }    
}

outdir <- normalizePath(outdir, mustWork=TRUE)

if (is.null(seg.file)) {
    tumor.coverage.file <- normalizePath(tumor.coverage.file, mustWork=TRUE)
}

if (is.null(sampleid)) stop("Need sampleid.")

library(PureCN)

PureCN(tumor.coverage.file, normal.coverage.file, tumor.vcf, genome,
gc.gene.file, seg.file=seg.file, snp.blacklist=snp.blacklist, 
stats.file=stats.file, target.weight.file=target.weight.file, 
normalDB=normalDB, sampleid=sampleid, post.optimize=post.optimize, 
normal.panel.vcf.file=normal.panel.vcf.file, pool=pool, 
model.homozygous=model.homozygous, outvcf=outvcf, 
outdir=outdir) 

