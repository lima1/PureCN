#' Power calculation for detecting somatic mutations
#'
#' This function calculates the probability of correctly rejecting the null
#' hypothesis that an alt allele is a sequencing error rather than a true
#' (mono-)clonal mutation.
#'
#'
#' @param coverage Mean sequencing coverage.
#' @param f Mean expected allelic fraction. If \code{NULL}, requires purity and
#' ploidy and then calculates the expected fraction.
#' @param purity Purity of sample. Only required when \code{f} is \code{NULL}.
#' @param ploidy Ploidy of sample. Only required when \code{f} is \code{NULL}.
#' @param cell.fraction Fraction of cells harboring mutation. Ignored if
#' \code{f} is not \code{NULL}.
#' @param error Estimated sequencing error rate.
#' @param fpr Required false positive rate for mutation vs. sequencing error.
#' @param verbose Verbose output.
#' @return A list with elements \item{power}{Power to detect somatic
#' mutations.} \item{k}{Minimum number of supporting reads.} \item{f}{Expected
#' allelic fraction. }
#' @author Markus Riester
#' @references Carter et al. (2012), Absolute quantification of somatic DNA
#' alterations in human cancer. Nature Biotechnology.
#' @examples
#'
#' purity <- c(0.1, 0.15, 0.2, 0.25, 0.4, 0.6, 1)
#' coverage <- seq(5, 35, 1)
#' power <- lapply(purity, function(p) sapply(coverage, function(cv)
#'     calculatePowerDetectSomatic(coverage = cv, purity = p, ploidy = 2,
#'     verbose = FALSE)$power))
#'
#' # Figure S7b in Carter et al.
#' plot(coverage, power[[1]], col = 1, xlab = "Sequence coverage",
#'     ylab = "Detection power", ylim = c(0, 1), type = "l")
#'
#' for (i in 2:length(power)) lines(coverage, power[[i]], col = i)
#' abline(h = 0.8, lty = 2, col = "grey")
#' legend("bottomright", legend = paste("Purity", purity),
#'     fill = seq_along(purity))
#'
#' # Figure S7c in Carter et al.
#' coverage <- seq(5, 350, 1)
#' power <- lapply(purity, function(p) sapply(coverage, function(cv)
#'     calculatePowerDetectSomatic(coverage = cv, purity = p, ploidy = 2,
#'         cell.fraction = 0.2, verbose = FALSE)$power))
#' plot(coverage, power[[1]], col = 1, xlab = "Sequence coverage",
#'     ylab = "Detection power", ylim = c(0, 1), type = "l")
#'
#' for (i in 2:length(power)) lines(coverage, power[[i]], col = i)
#' abline(h = 0.8, lty = 2, col = "grey")
#' legend("bottomright", legend = paste("Purity", purity),
#'     fill = seq_along(purity))
#'
#' @export calculatePowerDetectSomatic
calculatePowerDetectSomatic <- function(coverage, f = NULL,
purity = NULL, ploidy = NULL, cell.fraction = 1, error = 0.001, fpr = 5e-07,
verbose = TRUE) {
    # check parameters
    coverage <- round(coverage)
    if (coverage < 2) .stopUserError("coverage not in expected range (>=2)")
    if (error < 0 || error > 1) .stopUserError("error not in expected range.")
    if (fpr <= 0 || fpr > 1) .stopUserError("fpr not in expected range.")
    if (cell.fraction <= 0 || cell.fraction > 1) {
        .stopUserError("cell.fraction not in expected range.")
    }

    # calculate minimum number of required sequencing reads with mutation
    .pk <- function(m) {
        if (m == 0) return(1)
        dbinom(m, size = coverage, prob = error / 3)
     }
     k <- min(which(vapply(seq(0, coverage), .pk, double(1)) < fpr)) - 1
     if (verbose) message("Minimum ", k, " supporting reads.")

     # find allelic fraction to test
     if (is.null(f)) {
         if (is.null(purity) || is.null(ploidy)) {
             .stopUserError("Need either f or purity and ploidy.")
         }
         if (purity < 0 || purity > 1) {
            .stopUserError("purity not in expected range.")
         }
         if (ploidy <= 0) .stopUserError("ploidy not in expected range.")
         D <- purity * ploidy + (1 - purity) * 2
         f <- purity / D
         f <- f * cell.fraction
     }
     if (f < 0 || f > 1) .stopUserError("f not in expected range.")
     if (verbose) message("Expected allelic fraction ", f, ".")

     # calculate power
     .d <- function(k) {
         (fpr - .pk(k)) / (.pk(k - 1) - .pk(k))
     }
     power <- 1 - sum(dbinom(0:(k - 1), size = coverage, prob = f)) +
              .d(k) * dbinom(k, size = coverage, prob = f)

    list(
        power = power,
        k = k,
        f = f
    )
}
