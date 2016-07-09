createSNPBlacklist <- structure(function(# Create SNP black list
### Function to create a black list of germline SNPs with expected allelic
### fraction (AF) smaller than 0.5 in diploid genomes.  
vcf.files, 
### List of VCF files. When a VCF file contains multiple samples, 
### it will ignore all samples except the first.
n=min(10, length(vcf.files)),
### Required number of VCF files showing low allelic fraction to 
### blacklist a SNP id. 
low.af=0.025,
### Defines a low AF p-value.
high.af=0.1,
### Defines a high AF p-value. For every sample with high AF p-value, 
### there must be one more sample with low AF to reach the cutoff.  
chr.hash=NULL,
### Mapping of non-numerical chromsome names to numerical names
### (e.g. chr1 to 1, chr2 to 2, etc.). If NULL, assume chromsomes
### are properly ordered.   
genome="hg19"
### Version of the reference genome, required for the \code{readVcf} 
### function.
) {
    vcfs <- lapply(vcf.files, .readAndCheckVcf, genome)
    vcfs <- lapply(vcfs, function(x) x[info(x)$DB & 
        do.call(rbind, geno(x)$FA[,1, drop=FALSE])[,1]< 0.9 ,])

    .testAllelicRatioHeterozygousBias <- 
    function(vcf) {
        ar_all <- do.call(rbind, geno(vcf)$FA[,1, drop=FALSE])
        dp_all <- geno(vcf)$DP[,1, drop=FALSE]
        xx <-sapply(seq_along(ar_all), function(j) 
            pbeta(1/2, 
                shape1=ar_all[j]*dp_all[j]+1, 
                shape2=(1-ar_all[j])*dp_all[j]+1,log.p=FALSE))
    }
    vcfs.lh <- lapply(vcfs, .testAllelicRatioHeterozygousBias)
    vcfs.ar <- lapply(vcfs, function(vcf) do.call(rbind, geno(vcf)$FA[,1, drop=FALSE]))

    vcfs.smaller <- lapply(seq_along(vcfs), function(i) 
        vcfs[[i]][vcfs.lh[[i]] > 1 - low.af,])

    xx <- sort(table(do.call(c, lapply(vcfs.smaller, function(x) names(rowRanges(x))))))
    vcfs.greater <- lapply(seq_along(vcfs), function(i) 
        vcfs[[i]][vcfs.lh[[i]] < 1 - high.af,])
    xx.s <- sort(table(do.call(c, lapply(vcfs.greater, function(x) names(rowRanges(x))))))

    countTable <- data.frame(Count=as.vector(xx), row.names=names(xx))
    countTable <- cbind(countTable, Count.G=as.vector(xx.s[rownames(countTable)]))
    countTable$Count.G[is.na(countTable$Count.G)] <- 0

    snp.bl <- countTable[(countTable$Count-countTable$Count.G)>=n,,drop=FALSE]
    d.f <- do.call(rbind, lapply(vcfs, function(x) { 
        x <- x[rownames(x) %in% rownames(snp.bl)]
        data.frame(
            ID=rownames(x), 
            AR=do.call(rbind, geno(x)$FA[,1,drop=FALSE]),
            seqnames=as.character(seqnames(x)), 
            start=start(x)
        )
    }))

    mean.ar <- sapply(split(d.f$AR, d.f$ID), mean)
    snp.bl$Mean.AR <- mean.ar[rownames(snp.bl)]
    snp.bl$chr <- d.f$seqnames[match(rownames(snp.bl), d.f$ID)]
    snp.bl$start <- d.f$start[match(rownames(snp.bl), d.f$ID)]

    # segment
    d.f <- do.call(rbind, lapply(vcfs, function(x) 
        data.frame(binary=rownames(x) %in% rownames(snp.bl), 
                   seqnames=as.character(seqnames(x)), start=start(x))))

    d.f <- d.f[!duplicated(paste(d.f$seqnames, d.f$start)),]

    snp.bl.segmented <- segment(CNA(ifelse(d.f$binary,1,0), 
        chrom=d.f$seqnames, 
        maploc=d.f$start, 
        data.type="binary", 
        presorted=FALSE))$output

    if (is.null(chr.hash)) chr.hash <- .getChrHash(d.f$seqnames)

    snp.bl.segmented <- snp.bl.segmented[order(.strip.chr.name(snp.bl.segmented$chrom, chr.hash)),]
    
    ##value<< A list with elements
    list(
        snp.blacklist=snp.bl, ##<< A data.frame with blacklisted SNPs.
        segmented=snp.bl.segmented[,-1] ##<< A data.frame with blacklisted
## regions.
    )
##end<<
}, ex=function() {
# Assume VCF files of normals (for example obtained by a MuTect artifact
# detection run) are in directory poolofnormals:
mutect.normal.files <- dir("poolofnormals", pattern="vcf$", full.names=TRUE) 

# These files do not exist in our example, so we do not run the function here.
#snp.blacklist <- createSNPBlacklist(mutect.normal.files)
})
