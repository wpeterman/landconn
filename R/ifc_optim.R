# Optimize the incidence scale parameter (alpha) -------------------------------

#' Optimize the incidence function scale parameter (alpha)
#'
#' Estimate the dispersal scale `alpha` by finding the value whose incidence
#' function connectivity best predicts a response in a fitted model. This is a scale of
#' effect optimization: `alpha` is chosen to maximize the model's fit, where one
#' predictor is the connectivity produced by [ifc()]. Two ways to call it: pass a
#' `formula` and `data` to fit a GLM internally, or pass an already fitted model (a
#' `glm`, a negative binomial `glm.nb`, or an `unmarked` model) that contains a
#' `connectivity` term and have it refit across `alpha`.
#'
#' @param object Either a model `formula` (with `data` and `family`) or a fitted model to
#'   use as a template. As a formula, the predictor named `connectivity` is reserved and
#'   recomputed by [ifc()] at each candidate `alpha`, for example
#'   `occupied ~ connectivity + habitat`. As a fitted model (`glm`, `glm.nb`/`negbin`, or
#'   `unmarked`), it must already include a `connectivity` term (fit it once with a
#'   placeholder, such as `connectivity = 0`); [ifc_optim()] then refits it across `alpha`
#'   with [stats::update()]. For `unmarked`, `connectivity` lives in the `siteCovs`.
#' @param data A data frame with one row per site, in the same order as the rows of
#'   `dist_mat`. It holds the response and any covariates in the formula (but not
#'   `connectivity`, which is generated internally). Used only with the formula method.
#' @param dist_mat A square pairwise distance matrix, such as the output of [dist_mat()].
#'   It must have one row per site, in the same order as the model's sites.
#' @param model Numeric, one or more values in 1-4, the [ifc()] model(s) used to build
#'   connectivity. See [ifc()] for the four formulations. (Default = 1). A single value
#'   returns one fit; a vector (for example `1:4`) fits the best `alpha` for each model
#'   and returns a ranked model selection table. See Value and Details.
#' @param family A GLM family, as a family object, a family function, or its name (for
#'   example `binomial`, `poisson`, `gaussian`). (Default = `binomial()`). Used only with
#'   the formula method; a fitted template carries its own family.
#' @param ... Arguments passed between methods.
#' @param patch_area,occ_sites Optional numeric vectors passed through to [ifc()] (patch
#'   areas for Models 2-4, and a 0/1 occupancy vector to restrict contributors). Each
#'   must be equal in length to the number of sites. See [ifc()].
#' @param area_transform Passed to [ifc()]: `"none"` (Default) uses patch area directly,
#'   `"log"` uses `log(patch_area)`. You may supply both, `c("none", "log")`, to compare
#'   them for Models 2-4 in a model selection table (Model 1 ignores patch area, so it
#'   appears once). See [ifc()].
#' @param alpha_range Numeric length-2 vector giving the lower and upper bounds of the
#'   `alpha` search, in the same distance units as `dist_mat`. If `NULL` (Default), it is
#'   set from the spread of pairwise distances. If the estimate lands at a bound, widen
#'   this range and refit.
#' @param n_grid Number of log-spaced `alpha` values in the profile grid (Default = 50).
#'   The grid both seeds the optimizer and sets the resolution of the profile CI; raise
#'   it for a finer profile.
#' @param n_boot Number of parametric bootstrap replicates for the `alpha` confidence
#'   interval (Default = 0, no bootstrap). Each replicate simulates a new response from
#'   the fitted model, holding the site network fixed, and re-optimizes `alpha`. See
#'   Details.
#' @param conf_level Confidence level for the `alpha` interval(s) (Default = 0.95).
#'
#' @return When a single model and area transform are requested, an object of class
#'   `ifc_optim`. When several are requested (a vector `model`, or `area_transform =
#'   c("none", "log")` with Models 2-4), an object of class `ifc_optim_set`: a list with
#'   `fits` (the named list of individual `ifc_optim` fits), `table` (a ranked
#'   `ifc_modsel` model selection table, see [ifc_modsel()]), and `best` (the top-ranked
#'   `ifc_optim` fit). See [print.ifc_optim_set()].
#'
#'   A single `ifc_optim` object is a list with elements:
#'   \describe{
#'     \item{`alpha`}{The estimated scale parameter, in the distance units of `dist_mat`.
#'       Larger values mean connection probability decays more slowly with distance, so
#'       even distant sites stay connected.}
#'     \item{`alpha_ci`}{The profile likelihood confidence interval for `alpha` (the
#'       range of grid `alpha` values within `qchisq(conf_level, 1) / 2` log-likelihood
#'       units of the maximum). Its resolution is the grid; raise `n_grid` to sharpen it.}
#'     \item{`alpha_ci_boot`}{The bootstrap percentile interval for `alpha`, or `NULL`
#'       if `n_boot = 0`.}
#'     \item{`boot`}{The vector of bootstrap `alpha` estimates, or `NULL`.}
#'     \item{`logLik`, `AIC`}{The maximized log-likelihood and the AIC of the best model.
#'       `AIC` counts one extra parameter for the estimated `alpha`, so it is comparable
#'       across [ifc()] models and to a null (connectivity-free) model.}
#'     \item{`model`}{The fitted `glm` object at `alpha`. Treat its coefficients and
#'       standard errors as conditional on `alpha` (see Details).}
#'     \item{`profile`}{A data frame of `alpha` and `logLik` over the grid, used by
#'       `plot()`.}
#'     \item{`boundary`}{`TRUE` if `alpha` sits at a bound of `alpha_range`, a sign the
#'       range should be widened.}
#'   }
#'   Use [print.ifc_optim()], [summary.ifc_optim()], and [plot.ifc_optim()] to read it.
#'
#' @details
#'   `alpha` enters [ifc()] nonlinearly, so it cannot be read from a single GLM. For each
#'   candidate `alpha`, connectivity is recomputed, the GLM is fit, and its
#'   log-likelihood is recorded. The `alpha` that maximizes the log-likelihood is the
#'   estimate. This is a *profile likelihood*: a coarse log-spaced grid locates the peak
#'   and [stats::optimize()] refines it.
#'
#'   \strong{How to read the result.} Start with the profile plot. A sharp peak means
#'   `alpha` is well identified; a flat profile means the data carry little information
#'   about the scale, and the estimate should be treated with caution. The `alpha`
#'   estimate is a distance: it is the scale over which connection probability falls to
#'   `exp(-1)`, about 0.37, of its maximum.
#'
#'   \strong{Uncertainty.} Two intervals are available. The profile likelihood interval
#'   is fast and comes from the curvature of the profile. The parametric bootstrap
#'   (`n_boot > 0`) simulates new responses from the fitted model and re-optimizes
#'   `alpha` for each, giving a percentile interval. The bootstrap is parametric on
#'   purpose: connectivity is a whole-network quantity, so resampling sites would break
#'   the network (duplicated sites sit at zero distance). Holding the network fixed and
#'   resampling only the stochastic response is the coherent choice.
#'
#'   \strong{Caveats.} The GLM treats connectivity as a fixed predictor, so the GLM's own
#'   standard errors and p-values for the `connectivity` coefficient are conditional on
#'   the estimated `alpha` and are mildly anticonservative. If the response and
#'   `occ_sites` are the same occupancy data, you are predicting occupancy partly from
#'   itself; avoid that.
#'
#'   \strong{Model classes.} The formula method fits a GLM. To use another class, fit it
#'   once with a placeholder `connectivity` term and pass the fitted object: a `glm`, a
#'   negative binomial model from `MASS::glm.nb` (class `negbin`), or any `unmarked` model
#'   (for `unmarked`, put `connectivity` in the `siteCovs`). At each `alpha` the template
#'   is refit with [stats::update()], so its own parameters (the negative binomial theta,
#'   the detection parameters of an occupancy model) are re-estimated each step. Parameter
#'   counts in the information criteria include those parameters plus `alpha`. Keep the
#'   template simple: complex offsets or weights may not survive the refit.
#'
#'   \strong{Comparing models.} Passing several models (or both area transforms) fits the
#'   best `alpha` for each and ranks them in a model selection table. The table counts the
#'   number of parameters (`K`) correctly: the GLM coefficients, the dispersion parameter
#'   for families that estimate one (such as `gaussian`), and one more for the optimized
#'   `alpha`. Ranking is by AICc, with delta values and Akaike weights. The same parameter
#'   counting is used by the [AIC()], [BIC()], and [AICc()] methods for `ifc_optim`
#'   objects, so a fit compares correctly against an ordinary `glm` with no connectivity
#'   term. See [ifc_modsel()].
#'
#' @examples
#' set.seed(1)
#'
#' ## Simulate sites whose response depends on connectivity at a known alpha = 4
#' loc  <- data.frame(x = runif(120, 0, 20), y = runif(120, 0, 20))
#' d    <- dist_mat(loc)
#' conn <- as.numeric(ifc(4, d, model = 1, scale = TRUE))
#' dat  <- data.frame(abundance = 2 + 8 * conn + rnorm(120, 0, 0.5))
#'
#' ## Recover alpha (a continuous response identifies the scale well)
#' fit <- ifc_optim(abundance ~ connectivity, data = dat,
#'                  dist_mat = d, model = 1, family = gaussian)
#' fit
#' summary(fit)
#' plot(fit)
#'
#' ## Compare ifc models 1-4 (and both area transforms) in one call
#' p_a <- runif(120, 5, 50)
#' fits <- ifc_optim(abundance ~ connectivity, data = dat, dist_mat = d,
#'                   model = 1:4, patch_area = p_a,
#'                   area_transform = c("none", "log"), family = gaussian)
#' fits                 # ranked model selection table
#' fits$best            # the top-ranked single fit
#'
#' \donttest{
#' ## Add a parametric bootstrap confidence interval for alpha
#' fit_b <- ifc_optim(abundance ~ connectivity, data = dat,
#'                    dist_mat = d, model = 1, family = gaussian,
#'                    n_boot = 200)
#' fit_b$alpha_ci_boot
#'
#' ## Other model classes: fit a template with a placeholder `connectivity`, then pass it
#' if (requireNamespace("MASS", quietly = TRUE)) {
#'   nb_dat <- data.frame(counts = rnbinom(120, mu = exp(0.5 + 3 * conn), size = 3),
#'                        connectivity = 0)        # placeholder connectivity term
#'   nb0 <- MASS::glm.nb(counts ~ connectivity, data = nb_dat)
#'   ifc_optim(nb0, dist_mat = d, model = 1)
#' }
#' }
#'
#' @seealso [ifc()], [dist_mat()], [ifc_modsel()], [AICc()]
#'
#' @importFrom stats glm logLik optimize qchisq binomial simulate quantile AIC coef as.formula
#' @importFrom graphics abline rug
#' @export
#' @author Bill Peterman <peterman.73@@osu.edu>

ifc_optim <- function(object, ...) {
  UseMethod("ifc_optim")
}

#' @rdname ifc_optim
#' @exportS3Method ifc_optim formula
ifc_optim.formula <- function(object, data, dist_mat, model = 1,
                              family = stats::binomial(), patch_area = NULL,
                              occ_sites = NULL, area_transform = "none",
                              alpha_range = NULL, n_grid = 50, n_boot = 0,
                              conf_level = 0.95, ...) {

  formula <- stats::as.formula(object)
  if(!("connectivity" %in% all.vars(formula))){
    stop("`formula` must include the reserved predictor `connectivity`, e.g. `y ~ connectivity`.")
  }
  if(!is.data.frame(data)){
    stop("`data` must be a data frame with one row per site.")
  }

  ## Normalize the family argument (object, function, or name)
  if(is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
  if(is.function(family))  family <- family()

  adapter <- adapter_formula(formula, data, family)
  run_ifc_optim(adapter, dist_mat, model = model, patch_area = patch_area,
                occ_sites = occ_sites, area_transform = area_transform,
                alpha_range = alpha_range, n_grid = n_grid, n_boot = n_boot,
                conf_level = conf_level, call = match.call())
}


# Engine-agnostic core ---------------------------------------------------------

# An "adapter" is a list of closures that hides how a particular model class is fit
# and queried. Every adapter provides: n_sites, label, family, fit(conn, response),
# loglik(fitted), npar(fitted), simulate(fitted, nsim). The optimization below is
# written entirely against this interface, so it works for any engine.

# Adapter for the formula method: fit a GLM with connectivity as a predictor.
adapter_formula <- function(formula, data, family) {
  resp_name <- all.vars(formula[[2]])[1]
  list(
    n_sites = nrow(data),
    label   = paste0("glm (", family$family, ")"),
    family  = family$family,
    fit = function(conn, response = NULL) {
      dat <- data
      dat$connectivity <- conn
      if(!is.null(response)) dat[[resp_name]] <- response
      tryCatch(suppressWarnings(stats::glm(formula, data = dat, family = family)),
               error = function(e) NULL)
    },
    loglik   = function(fitted) as.numeric(stats::logLik(fitted)),
    npar     = function(fitted) as.integer(attr(stats::logLik(fitted), "df")),
    simulate = function(fitted, nsim) stats::simulate(fitted, nsim = nsim)
  )
}

# Optimize alpha for one (model, area_transform) combination using an adapter.
optimize_engine <- function(adapter, dist_mat, model, area_transform,
                            patch_area, occ_sites, alpha_range,
                            n_grid, n_boot, conf_level) {

  conn_at <- function(alpha) {
    as.numeric(ifc(alpha, dist_mat, model = model,
                   patch_area = patch_area, occ_sites = occ_sites,
                   area_transform = area_transform))
  }

  ll_at <- function(alpha, response = NULL) {
    fit <- adapter$fit(conn_at(alpha), response = response)
    if(is.null(fit)) return(NA_real_)
    adapter$loglik(fit)
  }

  ## grid scan, then refine around the best cell
  optim_alpha <- function(response = NULL) {
    grid <- exp(seq(log(alpha_range[1]), log(alpha_range[2]), length.out = n_grid))
    gll  <- vapply(grid, function(a) ll_at(a, response), numeric(1))
    if(all(is.na(gll))){
      stop("The model could not be fit at any alpha. Check the model, family, and data.")
    }
    k   <- which.max(gll)
    opt <- stats::optimize(function(a) -ll_at(a, response),
                           lower = grid[max(1, k - 1)], upper = grid[min(n_grid, k + 1)])
    list(alpha = opt$minimum, grid = grid, gll = gll, logLik = -opt$objective)
  }

  res  <- optim_alpha()
  best <- adapter$fit(conn_at(res$alpha))

  ## profile likelihood CI: grid alphas within the chi-square threshold of the max
  thresh   <- res$logLik - stats::qchisq(conf_level, 1) / 2
  in_ci    <- res$grid[!is.na(res$gll) & res$gll >= thresh]
  alpha_ci <- if(length(in_ci)) range(in_ci) else c(NA_real_, NA_real_)
  boundary <- res$alpha <= alpha_range[1] * 1.01 || res$alpha >= alpha_range[2] * 0.99

  df <- adapter$npar(best) + 1L          # model parameters plus alpha

  ## optional parametric bootstrap
  boot <- NULL
  alpha_ci_boot <- NULL
  if(n_boot > 0){
    sims <- tryCatch(adapter$simulate(best, n_boot), error = function(e) NULL)
    if(is.null(sims)){
      warning("Parametric bootstrap is unavailable for this model; returning the profile CI only.")
    } else {
      boot <- vapply(seq_len(n_boot), function(b) {
        tryCatch(optim_alpha(response = sims[[b]])$alpha, error = function(e) NA_real_)
      }, numeric(1))
      a_lo <- (1 - conf_level) / 2
      alpha_ci_boot <- stats::quantile(boot, c(a_lo, 1 - a_lo), na.rm = TRUE)
    }
  }

  structure(list(alpha = res$alpha,
                 alpha_ci = alpha_ci,
                 alpha_ci_boot = alpha_ci_boot,
                 boot = boot,
                 logLik = res$logLik,
                 df = df,
                 nobs = adapter$n_sites,
                 AIC = -2 * res$logLik + 2 * df,    # df already counts alpha
                 model = best,
                 profile = data.frame(alpha = res$grid, logLik = res$gll),
                 engine = adapter$label,
                 family = adapter$family,
                 ifc = list(model = model, area_transform = area_transform),
                 alpha_range = alpha_range,
                 conf_level = conf_level,
                 boundary = boundary,
                 call = NULL),
            class = "ifc_optim")
}

# Run the (possibly multi-model) optimization for any adapter, returning a single
# ifc_optim or, for several combinations, an ifc_optim_set.
run_ifc_optim <- function(adapter, dist_mat, model, patch_area, occ_sites,
                          area_transform, alpha_range, n_grid, n_boot,
                          conf_level, call) {

  if(!inherits(dist_mat, "matrix") || nrow(dist_mat) != ncol(dist_mat)){
    stop("`dist_mat` must be a square distance matrix.")
  }
  if(adapter$n_sites != nrow(dist_mat)){
    stop("The model must have the same number of sites as rows in `dist_mat`.")
  }
  model <- as.integer(model)
  if(any(is.na(model)) || any(!model %in% 1:4)){
    stop("`model` must be value(s) in 1:4. See `ifc()`.")
  }
  if(any(model > 1) && is.null(patch_area)){
    stop("`patch_area` is required for models 2, 3, and 4. See `ifc()`.")
  }
  area_transform <- match.arg(area_transform, c("none", "log"), several.ok = TRUE)

  ## default search range from the pairwise distances
  dvec <- lower(dist_mat)
  dpos <- dvec[dvec > 0]
  if(is.null(alpha_range)){
    alpha_range <- c(min(dpos) / 5, max(dvec) * 2)
  }
  if(length(alpha_range) != 2 || alpha_range[1] <= 0 || alpha_range[1] >= alpha_range[2]){
    stop("`alpha_range` must be a length-2 vector of positive, increasing values.")
  }

  ## build the grid of (model, area_transform) combinations
  ## Model 1 ignores patch area, so area_transform does not apply to it.
  combos <- do.call(rbind, lapply(unique(model), function(m) {
    if(m == 1L){
      data.frame(model = 1L, area_transform = "none", stringsAsFactors = FALSE)
    } else {
      data.frame(model = m, area_transform = unique(area_transform),
                 stringsAsFactors = FALSE)
    }
  }))
  combos <- unique(combos)

  fit_combo <- function(i) {
    f <- optimize_engine(adapter, dist_mat, combos$model[i], combos$area_transform[i],
                         patch_area, occ_sites, alpha_range, n_grid, n_boot, conf_level)
    f$call <- call
    f
  }

  ## single combination: return one ifc_optim object
  if(nrow(combos) == 1L){
    return(fit_combo(1))
  }

  ## multiple combinations: fit each and return an ifc_optim_set
  fits <- lapply(seq_len(nrow(combos)), fit_combo)
  names(fits) <- ifelse(combos$model == 1L, "model 1",
                        paste0("model ", combos$model, " (", combos$area_transform, ")"))
  tab  <- ifc_modsel(fits)
  best <- fits[[as.character(tab$model[1])]]

  structure(list(fits = fits, table = tab, best = best, call = call),
            class = "ifc_optim_set")
}
