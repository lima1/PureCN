test_that("test_createNormalDatabase", {
    normal.coverage.file <- system.file("extdata", "example_normal.txt", 
        package = "PureCN")
    normal2.coverage.file <- system.file("extdata", "example_normal2.txt", 
        package = "PureCN")
    normal3.coverage.file <- system.file("extdata", "example_normal3.cnn", 
        package = "PureCN")
    normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)
    normalDB <- createNormalDatabase(normal.coverage.files)
    expect_identical(normalDB$sex, c(NA, NA))
    expect_equal(findBestNormal(normal.coverage.files[1], normalDB), 
        normalizePath(normal.coverage.file[1]))
    pool <- findBestNormal(normal.coverage.files[1], normalDB, 
        num.normals = 2, pool = TRUE)
    n <- lapply(normal.coverage.files, readCoverageFile)
    expect_equal(length(pool), length(n[[1]]))
    expect_equal(normalDB$exon.median.coverage, apply(cbind(n[[1]]$average.coverage, 
        n[[2]]$average.coverage), 1, median))
    normalDB <- createNormalDatabase(normal.coverage.files, sex = c("A", 
        NA))
    expect_equal(normalDB$sex, as.character(c(NA, NA)))
    expect_equal(findBestNormal(normal.coverage.files[1], normalDB), 
        normalizePath(normal.coverage.file[1]))
    normalDB <- createNormalDatabase(normal.coverage.files, sex = c("A", 
        "F"))
    expect_equal(normalDB$sex, c(NA, "F"))
    expect_equal(normalDB$normal.coverage.files, 
                 sapply(normal.coverage.files, normalizePath, 
                        USE.NAMES = FALSE))

    expect_error(createNormalDatabase(normal.coverage.files, 
        sex = "A"))
    gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
        package = "PureCN")
    normal <- readCoverageFile(normal.coverage.file)
    correctCoverageBias(normal, gc.gene.file)
    suppressWarnings(correctCoverageBias(normal[sample(length(normal)), 
        ], gc.gene.file, "shuffled_gatk.txt"))
    createNormalDatabase(c(normal.coverage.files, "shuffled_gatk.txt"))
    tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
        package = "PureCN")
    best.normal.coverage.file <- findBestNormal(tumor.coverage.file, 
        normalDB)
    plotBestNormal(best.normal.coverage.file, tumor.coverage.file, 
        normalDB)
    expect_error(findBestNormal(normal3.coverage.file, normalDB, 
        num.normals = 2, pool = TRUE))
    expect_true(grepl("not align", geterrmessage()))
})

