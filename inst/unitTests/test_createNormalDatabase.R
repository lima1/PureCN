test_createNormalDatabase <- function() {
    gatk.normal.file <- system.file("extdata", "example_normal.txt", 
        package = "PureCN")
    gatk.normal2.file <- system.file("extdata", "example_normal2.txt", 
        package = "PureCN")
    gatk.normal.files <- c(gatk.normal.file, gatk.normal2.file)
    normalDB <- createNormalDatabase(gatk.normal.files)
    checkIdentical(c(NA, NA), normalDB$sex)
    checkEquals(normalizePath(gatk.normal.file[1]), 
        findBestNormal(gatk.normal.files[1], normalDB))

    normalDB <- createNormalDatabase(gatk.normal.files, sex=c("A", NA))
    checkEquals(as.character(c(NA, NA)), normalDB$sex)
    checkEquals(normalizePath(gatk.normal.file[1]), 
        findBestNormal(gatk.normal.files[1], normalDB))

    normalDB <- createNormalDatabase(gatk.normal.files, sex=c("A", "F"))
    checkEquals(c(NA, "F"), normalDB$sex)

    checkEquals(sapply(gatk.normal.files, normalizePath), 
        normalDB$gatk.normal.files, checkNames=FALSE)

    checkException(createNormalDatabase(gatk.normal.files, sex="A"), silent=TRUE) 
}    
