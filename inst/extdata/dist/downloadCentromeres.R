library(rtracklayer)
library(data.table)
library(PureCN)

data(chr.hash)
mySession <- browserSession("UCSC")
genomes <- c("hg18", "hg19", "hg38")
centromeres <- list()

for (genome in genomes) {
    genome(mySession) <- genome
    if (genome == "hg38") {
        tbl.gaps <- getTable( ucscTableQuery(mySession,track="Centromeres",
table="centromeres"))
    } else {
        tbl.gaps <- getTable( ucscTableQuery(mySession,  track="Gap",
            table="gap"))
        tbl.gaps <- tbl.gaps[tbl.gaps$type=="centromere",]
    }
    tbl.gaps.dt <- data.table(tbl.gaps)
    tbl.centromeres <- as.data.frame(tbl.gaps.dt[,
        list(chromStart=min(chromStart),chromEnd=max(chromEnd)),by=chrom])
    centromeres[[genome]] <- tbl.centromeres 
}

centromeres <- lapply(centromeres, function(x) {
    x$chromNumerical <- chr.hash$number[match(x$chrom, chr.hash$chr)]
    x[order(x$chromNumerical),1:3]
})

save(centromeres, file="data/centromeres.rda", compress="xz")
