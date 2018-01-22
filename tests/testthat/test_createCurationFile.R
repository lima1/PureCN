context("createCurationFile")

data(purecn.example.output)
file.rds <- tempfile(fileext = ".rds")
saveRDS(purecn.example.output, file = file.rds)

test_that("Example data is processed correctly", {
    ret <- createCurationFile(file.rds)
    expect_equal(ret$Purity, purecn.example.output$results[[1]]$purity)
    expect_equal(ret$Ploidy, purecn.example.output$results[[1]]$ploidy)
    expect_false(ret$Curated)
    expect_true(ret$Flagged)
    expect_equal(as.character(ret$Sampleid), purecn.example.output$input$sampleid)
})

test_that("Default curation file stores the first result", {
    retx <- readCurationFile(file.rds)
    expect_equal(retx$results[[1]]$purity, purecn.example.output$results[[1]]$purity)
    expect_equal(retx$results[[1]]$ploidy, purecn.example.output$results[[1]]$ploidy)
})

test_that("min.ploidy=2 ignores the first result", {
    retx <- readCurationFile(file.rds, min.ploidy = 2)
    expect_equal(retx$results[[1]]$purity, purecn.example.output$results[[2]]$purity)
    expect_equal(retx$results[[1]]$ploidy, purecn.example.output$results[[2]]$ploidy)
})

test_that("max.ploidy=2 ignores higher ploidy solutions", {
    retx <- readCurationFile(file.rds, max.ploidy = 2)
    expect_equal(sapply(retx$results, function(x) x$ploidy) < 
        2, rep(TRUE, length(retx$results)))
})

test_that("report.best.only works as expected", {
    retx <- readCurationFile(file.rds, report.best.only = TRUE)
    expect_equal(retx$results[[1]]$purity, purecn.example.output$results[[1]]$purity)
    expect_equal(retx$results[[1]]$ploidy, purecn.example.output$results[[1]]$ploidy)
    expect_equal(length(retx$results), 1)
})

test_that("overwriting works as expected", {
    retx <- purecn.example.output
    retx$results[[1]]$purity <- 0.8
    saveRDS(retx, file = file.rds)
    filename <- file.path(dirname(file.rds), paste0(gsub(".rds$", 
        "", basename(file.rds)), ".csv"))
    expect_warning(createCurationFile(file.rds, overwrite.uncurated = FALSE))
    ret <- read.csv(filename, as.is = TRUE)
    expect_equal(ret$Purity, purecn.example.output$results[[1]]$purity)
    expect_equal(ret$Ploidy, purecn.example.output$results[[1]]$ploidy)
    createCurationFile(file.rds)
    ret <- read.csv(filename, as.is = TRUE)
    expect_equal(ret$Purity, retx$results[[1]]$purity)
    expect_equal(ret$Ploidy, retx$results[[1]]$ploidy)
    ret$Curated <- TRUE
    write.csv(ret, file = filename, row.names = FALSE)
    saveRDS(purecn.example.output, file = file.rds)
    expect_warning(createCurationFile(file.rds))
    ret <- read.csv(filename, as.is = TRUE)
    expect_true(ret$Curated)
    expect_equal(ret$Purity, 0.8)
    ret$Ploidy <- 3.4
    write.csv(ret, file = filename, row.names = FALSE)
    retx <- readCurationFile(file.rds)
    expect_equal(ret$Purity, retx$results[[1]]$purity, tolerance=0.2)
    expect_equal(ret$Ploidy, retx$results[[1]]$ploidy, tolerance=0.5)
    ret$Purity <- "2.2w"
    write.csv(ret, file = filename, row.names = FALSE)
    expect_error(readCurationFile(file.rds))
    ret$Purity <- 2.2
    ret$Failed <- TRUE
    write.csv(ret, file = filename, row.names = FALSE)
    retx <- readCurationFile(file.rds, remove.failed = TRUE)
    expect_true(is.na(retx))
    ret$Failed <- "true"
    write.csv(ret, file = filename, row.names = FALSE)
    expect_error(readCurationFile(file.rds, remove.failed = TRUE), "logical")
    file.remove(filename)
})

test_that("warning occurs with missing curation file", {
    ret <- createCurationFile(file.rds)
    file.remove(gsub(".rds", ".csv", file.rds))
    expect_output(retx <- readCurationFile(file.rds), "does not exist, creating")
    expect_equal(retx$results[[1]]$purity, purecn.example.output$results[[1]]$purity)
    expect_equal(retx$results[[1]]$ploidy, purecn.example.output$results[[1]]$ploidy)
})

file.remove(file.rds)
