# S3 methods for ifc_optim_set objects ----------------------------------------

#' Methods for `ifc_optim_set` model selection objects
#'
#' An `ifc_optim_set` is returned by [ifc_optim()] when several [ifc()] models (or both
#' area transforms) are compared. It holds the individual fits, a ranked model selection
#' table, and the best fit.
#'
#' @param x,object An `ifc_optim_set` object.
#' @param ... Additional arguments passed to the underlying plotting calls (for `plot()`)
#'   or ignored.
#'
#' @return
#'   `print()` and `summary()` return their input invisibly after printing the model
#'   selection table; `summary()` also prints the best fit. `plot()` overlays the `alpha`
#'   log-likelihood profiles of all models and returns `NULL` invisibly.
#'
#' @details
#'   The best model is the top row of the table (lowest AICc by default). Drill into it
#'   with `object$best`, which is an ordinary [ifc_optim()] object with its own `print`,
#'   `summary`, and `plot` methods, or pull the table with `object$table`. See
#'   [ifc_modsel()] for the table columns.
#'
#' @examples
#' set.seed(1)
#' loc  <- data.frame(x = runif(120, 0, 20), y = runif(120, 0, 20))
#' d    <- dist_mat(loc)
#' conn <- as.numeric(ifc(4, d, model = 1, scale = TRUE))
#' dat  <- data.frame(abundance = 2 + 8 * conn + rnorm(120, 0, 0.5))
#' p_a  <- runif(120, 5, 50)
#'
#' fits <- ifc_optim(abundance ~ connectivity, data = dat, dist_mat = d,
#'                   model = 1:4, patch_area = p_a, family = gaussian)
#' fits
#' plot(fits)
#' fits$best
#'
#' @name ifc_optim_set-methods
#' @seealso [ifc_optim()], [ifc_modsel()]
NULL

#' @rdname ifc_optim_set-methods
#' @exportS3Method print ifc_optim_set
print.ifc_optim_set <- function(x, ...) {
  cat("Incidence scale model selection (`ifc_optim_set`)\n")
  cat("  ", length(x$fits), " models compared\n\n", sep = "")
  print(x$table)
  cat("\nBest model: ", x$table$model[1], "  (alpha = ",
      signif(x$best$alpha, 4), ")\n", sep = "")
  invisible(x)
}

#' @rdname ifc_optim_set-methods
#' @exportS3Method summary ifc_optim_set
summary.ifc_optim_set <- function(object, ...) {
  print(object)
  cat("\n--- Best model ---\n")
  print(object$best)
  invisible(object)
}

#' @rdname ifc_optim_set-methods
#' @exportS3Method plot ifc_optim_set
plot.ifc_optim_set <- function(x, ...) {
  fits <- x$fits
  finite_ll <- function(f) f$profile$logLik[is.finite(f$profile$logLik)]
  yr <- range(unlist(lapply(fits, finite_ll)))
  xr <- range(unlist(lapply(fits, function(f) f$profile$alpha)))
  cols <- seq_along(fits)

  plot(NA, xlim = xr, ylim = yr, log = "x",
       xlab = "alpha (dispersal scale)", ylab = "log-likelihood",
       main = "Incidence scale profiles by model", ...)
  for(i in seq_along(fits)){
    pr <- fits[[i]]$profile
    graphics::lines(pr$alpha, pr$logLik, col = cols[i], lwd = 2)
  }
  graphics::legend("bottomleft", legend = names(fits), col = cols,
                   lwd = 2, bty = "n", cex = 0.8)
  invisible(NULL)
}
