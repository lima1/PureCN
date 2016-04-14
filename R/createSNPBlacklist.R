createSNPBlacklist <- structure(function(# Create SNP black list
### Function to create a black list of germline SNPs with expected allelic
### fraction (AF) smaller than 0.5 in diploid genomes.  
vcf.files, 
### List of VCF files. When a VCF file contains multiple samples, it will
### ignore all samples except of the first.  
n=min(10, length(vcf.files)),
### Required number of VCF files showing low allelic fraction to blacklist a SNP id. 
low.af=0.025,
### Defines a low AF p-value.
high.af=0.1,
### Defines a high AF p-value. For every sample with high AF p-value, there
### must be one more sample with low AF to reach the cutoff.  
genome="hg19"
### Version of the reference genome, required for the readVcf() function.
) {
    vcfs <- lapply(vcf.files, readVcf, genome)
    vcfs <- lapply(vcfs, function(x) x[info(x)$DB & do.call(rbind, geno(x)$FA[,1, drop=FALSE])[,1]< 0.9 ,])

    .testLH <- 
    function(vcf) {
        ar_all <- do.call(rbind, geno(vcf)$FA[,1, drop=FALSE])
        dp_all <- geno(vcf)$DP[,1, drop=FALSE]
        xx <-sapply(1:length(ar_all), function(j) pbeta(1/2, shape1=ar_all[j]*dp_all[j]+1, shape2=(1-ar_all[j])*dp_all[j]+1,log.p=FALSE))
    }
    vcfs.lh <- lapply(vcfs, .testLH)
    vcfs.ar <- lapply(vcfs, function(vcf) do.call(rbind, geno(vcf)$FA[,1, drop=FALSE]))

    vcfs.smaller <- lapply(1:length(vcfs), function(i) vcfs[[i]][vcfs.lh[[i]] > 1 - low.af,])
    xx <- sort(table(do.call(c, lapply(vcfs.smaller, function(x) names(rowRanges(x))))))
    vcfs.greater <- lapply(1:length(vcfs), function(i) vcfs[[i]][vcfs.lh[[i]] < 1 - high.af,])
    xx.s <- sort(table(do.call(c, lapply(vcfs.greater, function(x) names(rowRanges(x))))))

    xx <- data.frame(Count=xx)
    xx <- cbind(xx, Count.G=xx.s[rownames(xx)])
    xx$Count.G[is.na(xx$Count.G)] <- 0

    snp.bl <- xx[(xx$Count-xx$Count.G)>=n,,drop=FALSE]
    d.f <- do.call(rbind, lapply(vcfs, function(x) { x <- x[rownames(x) %in% rownames(snp.bl)]; data.frame(ID=rownames(x), AR=do.call(rbind, geno(x)$FA[,1,drop=FALSE] ))}))
    mean.ar <- sapply(split(d.f$AR, d.f$ID), mean)
    snp.bl$Mean.AR <- mean.ar[rownames(snp.bl)]

    # segment
    d.f <- do.call(rbind, lapply(vcfs, function(x) data.frame(binary=rownames(x) %in% rownames(snp.bl), seqnames=as.character(seqnames(x)), start=start(x))))
    d.f <- d.f[!duplicated(paste(d.f$seqnames, d.f$start)),]

    snp.bl.segmented <- segment(CNA(ifelse(d.f$binary,1,0), chrom=d.f$seqnames, maploc=d.f$start, data.type="binary", presorted=FALSE))$output
    snp.bl.segmented <- snp.bl.segmented[order(.strip.chr.name(snp.bl.segmented$chrom)),]

    list(snp.black.list=snp.bl, segmented=snp.bl.segmented[,-1])
### A list with elements snp.black.list and segmented. 
### "snp.black.list" is just a list of SNP ids.
### "segmented" blacklists whole regions.
}, ex=function() {
# Assume VCF files of normals (for example obtained by a MuTect artifact
# detection run) are in directory poolofnormals:
mutect.normal.files <- dir("poolofnormals", pattern="vcf$", full.names=TRUE) 

# These files do not exist in our example, so we do not run the function here.
#snp.blacklist <- createSNPBlacklist(mutect.normal.files)
})
