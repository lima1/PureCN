test_calculatePowerDetectSomatic <- function() {
    p1 <- calculatePowerDetectSomatic(coverage=5, purity=1.0, ploidy=2)
    p2 <- calculatePowerDetectSomatic(coverage=5, f=0.5)
    checkEqualsNumeric(0.6407084, p1, tolerance=0.0001)
    checkEqualsNumeric(0.6407084, p2, tolerance=0.0001)
    p3 <- calculatePowerDetectSomatic(coverage=33, purity=0.5, ploidy=6)
    checkEqualsNumeric(0.8, p3, tolerance=0.001)
    p4 <- calculatePowerDetectSomatic(coverage=330, purity=0.2, ploidy=2, cell.fraction=0.2)
    checkEqualsNumeric(0.8, p4, tolerance=0.001)
}    

