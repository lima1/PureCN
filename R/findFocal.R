findFocal <- structure(function(# Find focal amplifications
### Function to find focal amplifications in segmented data. 
### This is automatically called in runAbsoluteCN.
seg,
### Segmentation data
size.cutoff=2000000,
### Cutoff for focal in base pairs.
cn.diff=2, 
### Minimum copy number delta between neighboring segments.
amp.cutoff=6
### Minimum amplification integer copy number.
) {
    focal <- rep(FALSE, nrow(seg))
    for (i in 1:nrow(seg)) {
        if (seg$C[i] < amp.cutoff) next
        if (seg$size[i] > size.cutoff) next
        size <- seg$size[i]
        if (i>1) {
            for (j in (i-1):1) {
                if (seg$C[j] < seg$C[i]-cn.diff) {
                    break
                }
                size <- size + seg$size[j]    
            }        
        }    
        if (i<nrow(seg)) {
            for (j in (i+1):nrow(seg)) {
                if (seg$C[j] < seg$C[i]-cn.diff) {
                    break
                }
                size <- size + seg$size[j]    
            }        
        }    
        focal[i] <- size < size.cutoff        
    }
    focal    
### Boolean vector for each segment wether it is focally 
### amplified or not.
},ex=function(){
gatk.normal.file <- system.file("extdata", "example_normal.txt", 
    package="PureCN")
gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
    package="PureCN")
vcf.file <- system.file("extdata", "example_vcf.vcf", 
    package="PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
    package="PureCN")

# Speed-up the runAbsoluteCN call by using the stored grid-search 
# (purecn.example.output$candidates).
data(purecn.example.output)

# The max.candidate.solutions parameter is set to a very low value only to
# speed-up this example.  This is not a good idea for real samples.
ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
    gatk.tumor.file=gatk.tumor.file,
    vcf.file=vcf.file, sampleid='Sample1', gc.gene.file=gc.gene.file,
    candidates=purecn.example.output$candidates, max.candidate.solutions=2,
    args.focal=list(size.cutoff = 2e+06), fun.focal=findFocal)
})
