test_createNormalDatabase <- function() {
    normal.coverage.file <- system.file("extdata", "example_normal.txt", 
        package = "PureCN")
    normal2.coverage.file <- system.file("extdata", "example_normal2.txt", 
        package = "PureCN")
    normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)
    normalDB <- createNormalDatabase(normal.coverage.files)
    checkIdentical(c(NA, NA), normalDB$sex)
    checkEquals(normalizePath(normal.coverage.file[1]), 
        findBestNormal(normal.coverage.files[1], normalDB))

    n <- lapply(normal.coverage.files, readCoverageGatk)

    checkEqualsNumeric(apply(cbind(n[[1]]$average.coverage, 
        n[[2]]$average.coverage),1,median), normalDB$exon.median.coverage)

    normalDB <- createNormalDatabase(normal.coverage.files, sex=c("A", NA))
    checkEquals(as.character(c(NA, NA)), normalDB$sex)
    checkEquals(normalizePath(normal.coverage.file[1]), 
        findBestNormal(normal.coverage.files[1], normalDB))

    normalDB <- createNormalDatabase(normal.coverage.files, sex=c("A", "F"))
    checkEquals(c(NA, "F"), normalDB$sex)

    checkEquals(sapply(normal.coverage.files, normalizePath), 
        normalDB$normal.coverage.files, checkNames=FALSE)

    checkException(createNormalDatabase(normal.coverage.files, sex="A"), silent=TRUE) 

    # create a GATK file with shuffled probes
    gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
        package = "PureCN")
    normal <- readCoverageGatk(normal.coverage.file)
    correctCoverageBias(normal, gc.gene.file)
    suppressWarnings(correctCoverageBias(normal[sample(nrow(normal)),], 
        gc.gene.file, "shuffled_gatk.txt"))
    checkException(createNormalDatabase(c(normal.coverage.files, 
        "shuffled_gatk.txt")))
    checkTrue(grepl("shuffled_gatk.txt", 
        geterrmessage()), msg=geterrmessage())
}    
