library('getopt')

### Parsing command line ------------------------------------------------------

spec <- matrix(c(
'help',           'h', 0, "logical",
'version',        'v', 0, "logical",
'force' ,         'f', 0, "logical",
'normal',         'n', 1, "character",
'tumor',          't', 1, "character",
'vcf',            'b', 1, "character",
'rds',            'r', 1, "character",
'genome',         'g', 1, "character",
'gcgene',         'c', 1, "character",
'segfile',        'l', 1, "character",
'snpblacklist',   's', 1, "character",
'normal_panel',   'p', 1, "character",
'statsfile',      'a', 1, "character",
'targetweightfile', 'e', 1, "character",
'normaldb',       'd', 1, "character",
'pool',           'm', 1, "integer",
'minaf',          'j', 1, "double",
'minpurity',      'k', 1, "double",
'maxpurity',      'x', 1, "double",
'postoptimize',   'z', 0, "logical",
'modelhomozygous','y', 0, "logical",
'outdir',         'o', 1, "character",
'outvcf',         'u', 0, "logical",
'sampleid',       'i', 1, "character"
), byrow=TRUE, ncol=4)
opt <- getopt(spec)

if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE))
    q(status=1)
}

if (!is.null(opt$version)) {
    message(as.character(packageVersion("PureCN")))
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
outvcf <- !is.null(opt$outvcf)
pool <- opt$pool
model.homozygous <- !is.null(opt$modelhomozygous)
file.rds <- opt$rds

if (!is.null(file.rds) && file.exists(file.rds)) {
    if (is.null(outdir)) outdir <- dirname(file.rds)
} else {
    if (is.null(sampleid)) stop("Need --sampleid.")
    if (is.null(genome)) stop("Need --genome")
    genome <- as.character(genome)
    file.rds <- file.path(outdir, paste0(sampleid, '_purecn.rds'))
    if (is.null(seg.file)) {
        tumor.coverage.file <- normalizePath(tumor.coverage.file, 
            mustWork=TRUE)
    }
}
    
outdir <- normalizePath(outdir, mustWork=TRUE)


library(PureCN)

### Run PureCN ----------------------------------------------------------------


if (file.exists(file.rds) && !force) {
    message(file.rds, 
        " already exists. Skipping... (--force will overwrite)")
    ret <- readCurationFile(file.rds)
    if (is.null(sampleid)) sampleid <- ret$input$sampleid
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
    file.log <- file.path(outdir, paste0(sampleid, '_purecn.log'))

    pdf(paste(outdir,"/", sampleid, '_purecn_segmentation.pdf', sep=''), 
        width=10, height=11)
    af.range = c(0.03, 0.97)
    if (!is.null(opt$minaf)) {
        af.range <- c(opt$minaf, 1-opt$minaf)
    }    
    test.purity <- seq(0.15, 0.95, by = 0.01)
    if (!is.null(opt$minpurity)) {
        if (!is.null(opt$maxpurity)) {
            test.purity <- seq(opt$minpurity, opt$maxpurity, by = 0.01)
        } else {
            test.purity <- seq(opt$minpurity, 0.95, by = 0.01)
        }
    }    
    ret <- runAbsoluteCN(normal.coverage.file=normal.coverage.file, 
            tumor.coverage.file=tumor.coverage.file, vcf.file=tumor.vcf,
            sampleid=sampleid, gc.gene.file=gc.gene.file, plot.cnv=TRUE,
            genome=genome, seg.file=seg.file,
            test.purity=test.purity,
            args.filterVcf=list(snp.blacklist=snp.blacklist, 
                af.range=af.range, stats.file=stats.file), 
            args.segmentation=list(target.weight.file=target.weight.file), 
            args.setMappingBiasVcf=
                list(normal.panel.vcf.file=normal.panel.vcf.file),
            normalDB=normalDB, model.homozygous=model.homozygous,
            model="betabin",
            log.file=file.log, post.optimize=post.optimize)
    dev.off()
    saveRDS(ret, file=file.rds)
}

### Create output files -------------------------------------------------------

createCurationFile(file.rds)
file.pdf <- file.path(outdir, paste0(sampleid, '_purecn.pdf'))
pdf(file.pdf, width=10, height=11)
plotAbs(ret, type='all')
dev.off()

file.png <- file.path(outdir, paste0(sampleid, '_purecn_contamination.png'))
png(file.png, width=800)
plotAbs(ret,1, type='contamination')
dev.off()

if (outvcf) {
    file.vcf <- file.path(outdir, paste0(sampleid, '_purecn.vcf'))
    vcfanno <- predictSomatic(ret, return.vcf=TRUE, 
        vcf.field.prefix="PureCN.")
    writeVcf(vcfanno, file=file.vcf)    
} else {
    file.csv <- file.path(outdir, paste0(sampleid, '_purecn_variants.csv'))
    write.csv(cbind(Sampleid=sampleid, predictSomatic(ret)), file=file.csv, 
        row.names=FALSE, quote=FALSE)
}    

file.loh <- file.path(outdir, paste0(sampleid, '_purecn_loh.csv'))
write.csv(cbind(Sampleid=sampleid, callLOH(ret)), file=file.loh, 
    row.names=FALSE, quote=FALSE)

file.genes <- file.path(outdir, paste0(sampleid, '_purecn_genes.csv'))
allAlterations <- callAlterations(ret, all.genes=TRUE)

write.csv(cbind(Sampleid=sampleid, gene.symbol=rownames(allAlterations), 
    allAlterations), row.names=FALSE, file=file.genes, quote=FALSE)

if (!is.null(ret$input$vcf)) {
    file.pdf <- file.path(outdir, paste0(sampleid, '_chromosomes_purecn.pdf'))
    pdf(file.pdf, width=9, height=10)
    vcf <- ret$input$vcf[ret$results[[1]]$SNV.posterior$beta.model$vcf.ids]
    chromosomes <- seqlevelsInUse(vcf)
    chromosomes <- chromosomes[orderSeqlevels(chromosomes)]
    for (chrom in chromosomes) {
        plotAbs(ret, 1, type='BAF', chr=chrom)
    }
    dev.off()
}
