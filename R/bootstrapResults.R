#' Bootstrapping variant fits
#'
#' This function bootstraps variants, then optionally re-ranks solutions by
#' using the bootstrap estimate of the likelihood score, and then optionally
#' removes solutions that never ranked high in any bootstrap replicate.
#'
#'
#' @param res Return object of the \code{\link{runAbsoluteCN}} function.
#' @param n Number of bootstrap replicates.
#' @param top Include solution if it appears in the top \code{n} solutions of
#' any bootstrap replicate. If \code{NULL}, do not filter solutions.
#' @param reorder Reorder results by bootstrap value.
#' @return Returns a \code{\link{runAbsoluteCN}} object with added bootstrap
#' value to each solution. This value
#' is the fraction of bootstrap replicates in which the solution ranked first.
#' @author Markus Riester
#' @seealso \code{\link{runAbsoluteCN}}
#' @examples
#'
#' data(purecn.example.output)
#' ret.boot <- bootstrapResults(purecn.example.output, n=100)
#' plotAbs(ret.boot, type="overview")
#'
#' @export bootstrapResults
#' @importFrom utils head
bootstrapResults <- function(res, n = 500, top = NULL, reorder = FALSE) {
    if (length(res$results) < 2) return(res)
    if (is.null(top)) top <- length(res$results)
    res$results <- .bootstrapResults(res$results, n = n, top = top,
        reorder = reorder)
    res
}

.bootstrapResults <- function(results, n, top, reorder) {
    ## Sample SNVs with replacement and recalculate log-likelihood.
    .bootstrapResult <- function(result) {
        lliks <- log(apply(result$SNV.posterior$likelihoods[
            !result$SNV.posterior$posteriors$FLAGGED, ], 1, max))
        lliks <- sum(sample(lliks, replace = TRUE))
        result$log.likelihood + sum(lliks) -
            sum(result$SNV.posterior$posteriors$FLAGGED)
    }
    best <- replicate(n, head(order(sapply(results, .bootstrapResult),
        decreasing = TRUE), top))

    ## Calculate bootstrap value as fraction solution is ranked first.
    bootstrap.value <- sapply(seq_along(results), function(i)
        sum(best[1, ] == i)) / ncol(best)
    for (i in seq_along(results)) {
        results[[i]]$bootstrap.value <- bootstrap.value[i]
    }

    ## Return only solutions that had ranked high in at least one replicate.
    best <- as.vector(best)
    results <- results[sort(unique(best))]
    if (reorder) {
        results <- results[order(sapply(results, function(x) x$bootstrap.value),
            decreasing = TRUE)]
    }
    .flagBootstrap(results)
}

.flagBootstrap <- function(results) {
    if (!is.null(results[[1]]$bootstrap.value)) {
        # max should be first, but be safe
        maxBootstrap <- max(sapply(results, function(r) r$bootstrap.value),
            na.rm = TRUE)
        if (maxBootstrap < 0.95) {
            for (i in seq_along(results)) {
                results[[i]]$flag <- TRUE
                results[[i]]$flag_comment <- .appendComment(
                    results[[i]]$flag_comment, "LOW BOOTSTRAP VALUE")
            }
        }
    }
    results
}
