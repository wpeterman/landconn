# Information criteria and model selection for ifc_optim -----------------------

# The point of these methods is to count parameters correctly. A fitted ifc_optim
# carries the GLM's own parameters (its coefficients, plus a dispersion parameter for
# families that estimate one, such as gaussian) and one additional parameter: the
# optimized `alpha`. `logLik.ifc_optim` reports that full degrees of freedom, so the
# base `AIC()`, `BIC()`, and the package's `AICc()` all count `alpha` automatically and
# compare correctly against an ordinary `glm`.

#' Log-likelihood and parameter count for an `ifc_optim` fit
#'
#' Returns the maximized log-likelihood with its degrees of freedom set to the GLM's
#' parameters plus one for the optimized `alpha`. This is what makes [AIC()], [BIC()],
#' and [AICc()] count the scale parameter.
#'
#' @param object An `ifc_optim` object from [ifc_optim()].
#' @param ... Ignored.
#' @return A `logLik` object. Its `df` attribute is the number of estimated parameters
#'   (GLM coefficients, any dispersion parameter, plus one for `alpha`) and its `nobs`
#'   attribute is the number of sites.
#' @seealso [AIC.ifc_optim()], [AICc()], [ifc_modsel()]
#' @importFrom stats AIC BIC logLik nobs
#' @exportS3Method logLik ifc_optim
logLik.ifc_optim <- function(object, ...) {
  ll  <- stats::logLik(object$model)
  val <- as.numeric(ll)
  attr(val, "df")   <- attr(ll, "df") + 1L      # GLM parameters plus alpha
  attr(val, "nobs") <- stats::nobs(object$model)
  class(val) <- "logLik"
  val
}

#' @rdname logLik.ifc_optim
#' @exportS3Method nobs ifc_optim
nobs.ifc_optim <- function(object, ...) stats::nobs(object$model)

#' Information criteria for `ifc_optim` fits
#'
#' `AIC()` and `BIC()` methods for [ifc_optim()] objects. They defer to the base
#' methods but, through [logLik.ifc_optim()], count one extra parameter for the
#' optimized `alpha`. Because of this you can compare an `ifc_optim` fit directly against
#' an ordinary `glm` (for example a null model with no connectivity term) with
#' `AIC(fit, null_glm)`.
#'
#' @param object An `ifc_optim` object.
#' @param ... Optional further fitted model objects (`ifc_optim` or `glm`) to compare, as
#'   in the base [stats::AIC()].
#' @param k The penalty per parameter, passed to [stats::AIC()] (Default = 2).
#' @return A numeric value, or a data frame of `df` and the criterion when several models
#'   are supplied, as for [stats::AIC()].
#' @seealso [AICc()], [ifc_modsel()], [logLik.ifc_optim()]
#' @exportS3Method AIC ifc_optim
AIC.ifc_optim <- function(object, ..., k = 2) NextMethod()

#' @rdname AIC.ifc_optim
#' @exportS3Method BIC ifc_optim
BIC.ifc_optim <- function(object, ...) NextMethod()


#' Second-order Akaike information criterion (AICc)
#'
#' Small-sample corrected AIC, \eqn{AICc = AIC + 2K(K + 1) / (n - K - 1)}, where `K` is
#' the number of parameters and `n` the number of observations. The default method works
#' for any model with a [logLik()] method, including [ifc_optim()] fits (for which `K`
#' includes the optimized `alpha`), `glm`, and `lm`.
#'
#' @param object A fitted model object with a `logLik` method.
#' @param ... Ignored.
#' @return The AICc value as a single number, or `NA` when `n - K - 1 <= 0` (too few
#'   observations for the correction to be defined).
#' @seealso [AIC.ifc_optim()], [ifc_modsel()]
#' @examples
#' set.seed(1)
#' loc  <- data.frame(x = runif(120, 0, 20), y = runif(120, 0, 20))
#' d    <- dist_mat(loc)
#' conn <- as.numeric(ifc(4, d, model = 1, scale = TRUE))
#' dat  <- data.frame(abundance = 2 + 8 * conn + rnorm(120, 0, 0.5))
#' fit  <- ifc_optim(abundance ~ connectivity, data = dat, dist_mat = d,
#'                   model = 1, family = gaussian)
#' AICc(fit)
#' @export
AICc <- function(object, ...) UseMethod("AICc")

#' @rdname AICc
#' @exportS3Method AICc default
AICc.default <- function(object, ...) {
  ll <- stats::logLik(object)
  k  <- attr(ll, "df")
  n  <- attr(ll, "nobs")
  if(is.null(n)) n <- stats::nobs(object)
  aic <- -2 * as.numeric(ll) + 2 * k
  if(n - k - 1 <= 0) return(NA_real_)
  aic + (2 * k * (k + 1)) / (n - k - 1)
}


#' Model selection table for `ifc_optim` fits
#'
#' Rank two or more [ifc_optim()] fits (or the fits inside an `ifc_optim_set`) by an
#' information criterion, with delta values and Akaike weights. The number of parameters
#' (`K`) counts the GLM coefficients, any dispersion parameter, and the optimized `alpha`.
#'
#' @param ... Two or more `ifc_optim` objects, a single named list of them, or an
#'   `ifc_optim_set`.
#' @param rank The information criterion to rank by: `"AICc"` (Default), `"AIC"`, or
#'   `"BIC"`. `delta` and `weight` are computed from the chosen criterion.
#' @return A data frame of class `ifc_modsel`, one row per model, sorted best first, with
#'   columns: `model` (label), `K` (parameter count), `alpha`, `logLik`, `AIC`, `AICc`,
#'   `BIC`, `delta` (difference from the best model on the ranking criterion), and
#'   `weight` (Akaike weight). Models within about 2 `delta` units of the best have
#'   comparable support.
#' @details Comparisons assume the fits share the same response and data. The function
#'   warns if the number of observations differs across fits. Weights are
#'   \eqn{w_i = \exp(-\Delta_i / 2) / \sum_j \exp(-\Delta_j / 2)} and sum to 1.
#' @seealso [ifc_optim()], [AICc()], [AIC.ifc_optim()]
#' @examples
#' set.seed(1)
#' loc  <- data.frame(x = runif(120, 0, 20), y = runif(120, 0, 20))
#' d    <- dist_mat(loc)
#' conn <- as.numeric(ifc(4, d, model = 1, scale = TRUE))
#' dat  <- data.frame(abundance = 2 + 8 * conn + rnorm(120, 0, 0.5))
#' p_a  <- runif(120, 5, 50)
#'
#' m1 <- ifc_optim(abundance ~ connectivity, data = dat, dist_mat = d, model = 1,
#'                 family = gaussian)
#' m2 <- ifc_optim(abundance ~ connectivity, data = dat, dist_mat = d, model = 2,
#'                 patch_area = p_a, family = gaussian)
#' ifc_modsel(m1, m2)
#' @export
ifc_modsel <- function(..., rank = c("AICc", "AIC", "BIC")) {
  rank <- match.arg(rank)
  dots <- list(...)

  ## Accept an ifc_optim_set, a single list of fits, or loose fits in `...`
  if(length(dots) == 1 && inherits(dots[[1]], "ifc_optim_set")){
    fits <- dots[[1]]$fits
  } else if(length(dots) == 1 && is.list(dots[[1]]) && !inherits(dots[[1]], "ifc_optim")){
    fits <- dots[[1]]
  } else {
    fits <- dots
  }
  if(!length(fits) || !all(vapply(fits, inherits, logical(1), "ifc_optim"))){
    stop("`ifc_modsel()` expects `ifc_optim` objects (or a list of them).")
  }

  ## Labels: use list names, else build from the ifc model and transform
  labs <- names(fits)
  if(is.null(labs)) labs <- rep("", length(fits))
  for(i in seq_along(fits)){
    if(!nzchar(labs[i])){
      f <- fits[[i]]
      labs[i] <- if(f$ifc$model == 1) "model 1" else
        paste0("model ", f$ifc$model, " (", f$ifc$area_transform, ")")
    }
  }

  n <- vapply(fits, stats::nobs, numeric(1))
  if(length(unique(n)) > 1){
    warning("Fits have different numbers of observations; the comparison may be invalid.")
  }

  K     <- vapply(fits, function(f) attr(stats::logLik(f), "df"), numeric(1))
  llv   <- vapply(fits, function(f) as.numeric(stats::logLik(f)), numeric(1))
  aic   <- vapply(fits, stats::AIC, numeric(1))
  aicc  <- vapply(fits, AICc, numeric(1))
  bic   <- vapply(fits, stats::BIC, numeric(1))
  alpha <- vapply(fits, function(f) f$alpha, numeric(1))

  ic    <- switch(rank, AICc = aicc, AIC = aic, BIC = bic)
  delta <- ic - min(ic, na.rm = TRUE)
  w     <- exp(-0.5 * delta) / sum(exp(-0.5 * delta), na.rm = TRUE)

  tab <- data.frame(model  = labs,
                    K      = K,
                    alpha  = round(alpha, 3),
                    logLik = round(llv, 2),
                    AIC    = round(aic, 2),
                    AICc   = round(aicc, 2),
                    BIC    = round(bic, 2),
                    delta  = round(delta, 2),
                    weight = round(w, 3),
                    stringsAsFactors = FALSE)
  tab <- tab[order(ic), ]
  rownames(tab) <- NULL
  structure(tab, class = c("ifc_modsel", "data.frame"), rank = rank)
}

#' @rdname ifc_modsel
#' @param x An `ifc_modsel` table.
#' @exportS3Method print ifc_modsel
print.ifc_modsel <- function(x, ...) {
  cat("Incidence scale model selection (ranked by ", attr(x, "rank"), ")\n", sep = "")
  cat("K counts GLM parameters plus the optimized alpha.\n\n")
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}
