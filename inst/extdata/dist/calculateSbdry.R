library(PureCN)

alpha <- formals(segmentationCBS)$alpha
eta <- formals(segment)$eta
nperm <- formals(segment)$nperm
max.ones <- floor(nperm * alpha) + 1
set.seed(123)

purecn.DNAcopy.bdry <- getbdry(eta, nperm, max.ones)
save(purecn.DNAcopy.bdry, file="~/git/PureCN/data/purecn.DNAcopy.bdry.rda", compress="xz")
