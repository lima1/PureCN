poolCoverage <- structure(function(
### Averages the coverage of a list of samples.
all.data, 
### List of normals, read with readCoverageGatk.
remove.chrs=c(),
### Remove these chromosomes from the pool.
w=NULL
### List of weights for each sample
) { 
    pool = all.data[[1]]
    if (length(all.data) == 1) {
        return(.removeChr(pool, remove.chrs))
    }
    if (is.null(w)) w <- rep(1,length(all.data))
    w <- w/w[1]

    for (i in 2:length(all.data)) {
#        pool$sequenced.base = find.max.of.2lists(pool$sequenced.base, 
#            all.data[[i]]$sequenced.base)
        pool$coverage <- pool$coverage + (w[i] * all.data[[i]]$coverage)
        pool$average.coverage <- pool$average.coverage + (w[i] * all.data[[i]]$average.coverage)
        pool$base.with..10.coverage <- pool$base.with..10.coverage + (w[i] * all.data[[i]]$base.with..10.coverage)
    }
    return(.removeChr(pool, remove.chrs))
### A data.frame with the averaged coverage over all normals.
},ex=function() {
gatk.normal.file <- system.file("extdata", "example_normal.txt", package="PureCN")
gatk.normal2.file <- system.file("extdata", "example_normal2.txt", package="PureCN")
gatk.normal.files <- c(gatk.normal.file, gatk.normal2.file)
gatk.tumor.file <- system.file("extdata", "example_tumor.txt", package="PureCN")

normalDB <- createNormalDatabase(gatk.normal.files)

# get the best 2 normals and average them
gatk.best.normal.files <- findBestNormal(gatk.tumor.file, normalDB, 
    num.normals=2)
pool <- poolCoverage(lapply(gatk.best.normal.files, readCoverageGatk),
     remove.chrs=c('chrX', 'chrY'))
})    

.removeChr <- function(pool, remove.chrs=c()) {
    pool$chr <- gsub("^chrchr", "chr", pool$chr)
    idx <- pool$chr %in% remove.chrs
    pool$coverage[idx] <- NA
    pool$average.coverage[idx] <- NA
    pool$base.with..10.coverage[idx] <- NA
    pool
}
