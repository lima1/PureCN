context("createNormalDatabase")

tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
    package = "PureCN")
normal.coverage.file <- system.file("extdata", "example_normal.txt", 
    package = "PureCN")
normal2.coverage.file <- system.file("extdata", "example_normal2.txt", 
    package = "PureCN")
normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)
normalDB <- createNormalDatabase(normal.coverage.files)

test_that("NormalDB of example data matches expectated values", {
    expect_identical(normalDB$sex, c(NA, NA))
    pool <- calculateTangentNormal(normal.coverage.files[1], normalDB)

    n <- lapply(normal.coverage.files, readCoverageFile)
    expect_equal(length(pool), length(n[[1]]))
    expect_equal(as.character(n[[1]]), normalDB$intervals)
 })

test_that("Provided sex is handled correctly", {
    expect_warning(
        normalDB2 <- createNormalDatabase(normal.coverage.files, sex = c("A",
            NA))
    )
    expect_equal(normalDB2$sex, as.character(c(NA, NA)))
    expect_warning(
        normalDB2 <- createNormalDatabase(normal.coverage.files, sex = c("A",
        "F"))
    )
    expect_equal(normalDB2$sex, c(NA, "F"))
    expect_equal(normalDB2$normal.coverage.files, 
                 sapply(normal.coverage.files, normalizePath, 
                        USE.NAMES = FALSE))

    expect_error(createNormalDatabase(normal.coverage.files, sex = "A"))
})

test_that("Exceptions happen with wrong input", {
    interval.file <- system.file("extdata", "example_intervals.txt", 
        package = "PureCN")
    normal <- readCoverageFile(normal.coverage.file)
    correctCoverageBias(normal, interval.file)
    output.file <- tempfile(fileext = ".txt")
    expect_output(correctCoverageBias(normal[sample(length(normal)), 
        ], interval.file, output.file), "WARN")
    createNormalDatabase(c(normal.coverage.files, output.file))
    best.normal.coverage.file <- calculateTangentNormal(tumor.coverage.file, 
        normalDB)
    normal3.coverage.file <- system.file("extdata", "example_normal3.cnn", 
        package = "PureCN")
    expect_error(calculateTangentNormal(normal3.coverage.file, normalDB),
       "not align")
    expect_error(createNormalDatabase(normal.coverage.file), "At least 2")
    expect_output(createNormalDatabase( c(normal.coverage.file, normal.coverage.file, 
                                          normal2.coverage.file)), "duplicated")
    file.remove(output.file)
})


test_that("Exceptions happen with outdated databases", {
    normalDB2 <- normalDB
    normalDB2$version <- NULL
    expect_error(runAbsoluteCN(normal.coverage.file, tumor.coverage.file, normalDB = normalDB2),
                 "incompatible")
    expect_error( calculateTangentNormal(tumor.coverage.file, normalDB2), "incompatible")
})


test_that("Exception thrown when user mixed gc-normalized and raw coverages.", {
    normal.coverage.files.wrong <- c(tempfile(fileext="_coverage.txt"), tempfile(fileext="_loess.txt"))
    file.create(normal.coverage.files.wrong)
    expect_error( createNormalDatabase(normal.coverage.files.wrong), "suffix")
    file.remove(normal.coverage.files.wrong)
})
