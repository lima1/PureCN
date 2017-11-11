context("calculatePowerDetectSomatic")

test_that("Power is calculated correctly for examples", {
    p1 <- calculatePowerDetectSomatic(coverage = 5, purity = 1, 
        ploidy = 2)$power
    p2 <- calculatePowerDetectSomatic(coverage = 5, f = 0.5)$power
    expect_equal(p1, 0.6407084, tolerance=0.0001)
    expect_equal(p2, 0.6407084, tolerance=0.0001)
    p3 <- calculatePowerDetectSomatic(coverage = 33, purity = 0.5, 
        ploidy = 6)$power
    expect_equal(p3, 0.8, tolerance=0.001)
    p4 <- calculatePowerDetectSomatic(coverage = 330, purity = 0.2, 
        ploidy = 2, cell.fraction = 0.2)$power
    expect_equal(p4, 0.8, tolerance=0.001)
})

test_that("Exceptions happen with wrong input", {
    expect_error(calculatePowerDetectSomatic(coverage = 5))
    expect_error(calculatePowerDetectSomatic(coverage = 5, f = 1.1))
    expect_error(calculatePowerDetectSomatic(coverage = 1, f = 0.9))
    expect_error(calculatePowerDetectSomatic(coverage = 3, purity = 1.1, 
        ploidy = 2))
    expect_error(calculatePowerDetectSomatic(coverage = 3, purity = 1, 
        ploidy = -1))
    expect_error(calculatePowerDetectSomatic(coverage = 5, cell.fraction = 1.1))
})
