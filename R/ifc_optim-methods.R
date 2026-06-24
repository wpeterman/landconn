# S3 methods for ifc_optim objects --------------------------------------------

#' Methods for `ifc_optim` objects
#'
#' Display, summarize, and plot the result of [ifc_optim()].
#'
#' @param x,object An object of class `ifc_optim`, as returned by [ifc_optim()].
#' @param ... Additional arguments passed to the underlying plotting calls (for
#'   `plot()`) or ignored (for `print()` and `summary()`).
#'
#' @return
#'   `print()` and `summary()` return their input invisibly after printing.
#'   `summary()` additionally prints the fitted GLM coefficient table. `plot()` draws the
#'   `alpha` profile (log-likelihood against `alpha`) with the estimate, the confidence
#'   threshold, and the interval marked, and returns `NULL` invisibly.
#'
#' @details
#'   The profile plot is the main diagnostic. A sharp peak means `alpha` is well
#'   identified by the data; a flat profile means it is not, and the estimate should be
#'   treated with caution. The dashed horizontal line is the confidence threshold, and
#'   the points where the profile crosses it bound the profile confidence interval.
#'
#' @examples
#' set.seed(1)
#' loc  <- data.frame(x = runif(120, 0, 20), y = runif(120, 0, 20))
#' d    <- dist_mat(loc)
#' conn <- as.numeric(ifc(4, d, model = 1, scale = TRUE))
#' dat  <- data.frame(abundance = 2 + 8 * conn + rnorm(120, 0, 0.5))
#' fit  <- ifc_optim(abundance ~ connectivity, data = dat, dist_mat = d,
#'                   model = 1, family = gaussian)
#'
#' fit
#' summary(fit)
#' plot(fit)
#'
#' @name ifc_optim-methods
#' @seealso [ifc_optim()]
NULL

#' @rdname ifc_optim-methods
#' @exportS3Method print ifc_optim
print.ifc_optim <- function(x, ...) {
  cat("Incidence scale (alpha) optimization (`ifc_optim`)\n")
  cat("  ifc model:   ", x$ifc$model,
      " | family: ", x$family, "\n", sep = "")
  cat("  alpha:       ", signif(x$alpha, 4), "\n", sep = "")
  cat("  ", round(100 * x$conf_level), "% profile CI: ",
      signif(x$alpha_ci[1], 4), " to ", signif(x$alpha_ci[2], 4), "\n", sep = "")
  if(!is.null(x$alpha_ci_boot)){
    cat("  ", round(100 * x$conf_level), "% bootstrap CI: ",
        signif(x$alpha_ci_boot[1], 4), " to ", signif(x$alpha_ci_boot[2], 4),
        " (", length(x$boot), " reps)\n", sep = "")
  }
  cat("  logLik:      ", round(x$logLik, 3),
      " | AIC: ", round(x$AIC, 2), " (includes alpha)\n", sep = "")
  if(isTRUE(x$boundary)){
    cat("\n  NOTE: alpha is at a bound of the search range; widen `alpha_range` and refit.\n")
  }
  invisible(x)
}

#' @rdname ifc_optim-methods
#' @exportS3Method summary ifc_optim
summary.ifc_optim <- function(object, ...) {
  print(object)
  cat("\nFitted GLM at the optimized alpha:\n")
  print(stats::coef(summary(object$model)))
  cat("\n(Standard errors above are conditional on the estimated alpha; see ?ifc_optim.)\n")
  invisible(object)
}

#' @rdname ifc_optim-methods
#' @exportS3Method plot ifc_optim
plot.ifc_optim <- function(x, ...) {
  prof <- x$profile[is.finite(x$profile$logLik), ]
  plot(prof$alpha, prof$logLik, type = "l",
       log = "x",
       xlab = "alpha (dispersal scale)",
       ylab = "log-likelihood",
       main = "Incidence scale profile",
       ...)
  ## Estimate and profile confidence threshold
  graphics::abline(v = x$alpha, col = "red", lwd = 2)
  thresh <- x$logLik - stats::qchisq(x$conf_level, 1) / 2
  graphics::abline(h = thresh, lty = 2, col = "grey40")
  if(all(is.finite(x$alpha_ci))){
    graphics::abline(v = x$alpha_ci, lty = 3, col = "red")
  }
  ## Bootstrap estimates, if present, as a rug along the x axis
  if(!is.null(x$boot)){
    graphics::rug(x$boot[is.finite(x$boot)], col = "steelblue")
  }
  invisible(NULL)
}
