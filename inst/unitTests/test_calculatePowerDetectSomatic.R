test_calculatePowerDetectSomatic <- function() {
    p1 <- calculatePowerDetectSomatic(coverage=5, purity=1.0, ploidy=2)$power
    p2 <- calculatePowerDetectSomatic(coverage=5, f=0.5)$power
    checkEqualsNumeric(0.6407084, p1, tolerance=0.0001)
    checkEqualsNumeric(0.6407084, p2, tolerance=0.0001)
    p3 <- calculatePowerDetectSomatic(coverage=33, purity=0.5, ploidy=6)$power
    checkEqualsNumeric(0.8, p3, tolerance=0.001)
    p4 <- calculatePowerDetectSomatic(coverage=330, purity=0.2, ploidy=2, cell.fraction=0.2)$power
    checkEqualsNumeric(0.8, p4, tolerance=0.001)
    # test wrong parameters
    checkException(calculatePowerDetectSomatic(coverage=5))
    checkException(calculatePowerDetectSomatic(coverage=5, f=1.1))
    checkException(calculatePowerDetectSomatic(coverage=1, f=0.9))
    checkException(calculatePowerDetectSomatic(coverage=3, purity=1.1, 
        ploidy=2))
    checkException(calculatePowerDetectSomatic(coverage=3, purity=1, 
        ploidy=-1))
    checkException(calculatePowerDetectSomatic(coverage=5, cell.fraction=1.1))
}    

