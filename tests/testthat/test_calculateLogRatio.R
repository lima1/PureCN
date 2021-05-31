context("calculateLogRatio")

test_that("Misaligned on- and off-target regions are aligned", {
    x <- readCoverageFile(
        system.file("extdata", "example_intervals_tiny_ot.txt.gz",
        package = "PureCN"))
    set.seed(123)
    l1 <- rnorm(length(x), mean = 0.25, sd=0.3)
    l2 <- rnorm(length(x), mean = -0.25, sd=0.3)
    x$log.ratio <- l1
    x$log.ratio[x$on.target] <- l2[x$on.target]
    expect_lt(t.test( x$log.ratio[x$on.target], x$log.ratio[!x$on.target])$p.value, 0.001)

    xc <- x
    xc$log.ratio <- PureCN:::.calibrate_off_target_log_ratio(x)
    expect_gt(t.test( xc$log.ratio[x$on.target], xc$log.ratio[!x$on.target])$p.value, 0.001)

    x$log.ratio <- l2
    x$log.ratio[x$on.target] <- l1[x$on.target]
    expect_lt(t.test( x$log.ratio[x$on.target], x$log.ratio[!x$on.target])$p.value, 0.001)

    xc <- x
    xc$log.ratio <- PureCN:::.calibrate_off_target_log_ratio(x)
    expect_gt(t.test( xc$log.ratio[x$on.target], xc$log.ratio[!x$on.target])$p.value, 0.001)

})

