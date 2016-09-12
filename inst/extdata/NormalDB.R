library('getopt')

spec <- matrix(c(
'help' , 'h', 0, "logical",
'intervalfile', 'i', 1, "character",
'referencefile', '', 1, "character",
'coveragefiles', 'c', 1, "character",
'vcffiles', 'v', 1, "character",
'outdir' , 'o', 1, "character"
), byrow=TRUE, ncol=4)
opt <- getopt(spec)

if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE))
    q(status=1)
}

if ( !is.null(opt$intervalfile ) {
    if (is.null(opt$referencefile)) {
        stop("Need reference Fasta file for creating gc.gene.file.")
    }
    calculateGCContentByInterval(opt$intervalfile, opt$referencefile,
        file.path(outdir,"gc_gene.txt"))
}

.checkFileList <- function(file) {
    files <- read.delim(file, as.is=TRUE, header=FALSE)[,1]
    numExists <- sum(file.exists(files), na.rm=TRUE)
    if (numExists < length(files)) { 
        stop("File not exists in file ", file)
    }
    files
}

if (is.null(opt$coveragefiles)) {
    stop("need --coveragefiles.")
}

coverageFiles <- .checkFileList(opt$coveragefiles)
outdir <- opt$outdir
if (is.null(outdir)) {
    stop("need --outdir")
}
outdir <- normalizePath(outdir, mustWork=TRUE)

if (length(coverageFiles)) {
    library(PureCN)
    normalDB <- createNormalDatabase(coverageFiles)
    saveRDS(normalDB, file=file.path(outdir,"normalDB.rds"))
}

if (length(coverageFiles) > 3) {
    library(PureCN)
    target.weight.file <- file.path(outdir,"target_weights.txt")
    createTargetWeights(coverageFiles[1:2], coverageFiles[-(1:2)], 
        target.weight.file)
} else {
    message("Not enough coverage files for creating target_weights.txt")
}


if (!is.null(opt$vcffiles)) {
    vcfFiles <- .checkFileList(opt$vcffiles)

    snp.blacklist <- file.path(outdir, c("SNP_blacklist.csv",
        "SNP_blacklist_segmented.csv"))

    if (length(vcfFiles) > 3) {
        library(PureCN)
        snp.bl <- createSNPBlacklist(vcfFiles)
        write.csv(snp.bl[[1]], file=snp.blacklist[1])
        write.csv(snp.bl[[2]], file=snp.blacklist[2], row.names=FALSE, 
           quote=FALSE)
    } else {
        message("Not enough VCF files for creating SNP_blacklist.csv")
    }
}
