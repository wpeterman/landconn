# Optimize the incidence scale parameter (alpha) -------------------------------

#' Optimize the incidence function scale parameter (alpha)
#'
#' Estimate the association scale `alpha` by finding the value whose incidence
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
#'   supported `unmarked`), it must already include a `connectivity` term (fit it once with a
#'   placeholder, such as `connectivity = 0`); [ifc_optim()] then refits it across `alpha`
#'   while preserving its fitted call settings. For supported single-season
#'   `unmarked::occu()` and `unmarked::pcount()` models, `connectivity` must occur
#'   in the state formula and in `siteCovs`, never in the detection formula.
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
#' @param n_grid Integer number of log-spaced `alpha` values in the profile grid
#'   (default 50, minimum 3).
#'   The grid both seeds the optimizer and sets the resolution of the profile CI; raise
#'   it for a finer profile.
#' @param n_boot Nonnegative integer number of parametric bootstrap replicates for the `alpha` confidence
#'   interval (Default = 0, no bootstrap). Each replicate simulates a new response from
#'   the fitted model, holding the site network fixed, and re-optimizes `alpha`. See
#'   Details.
#' @param conf_level Numeric confidence level strictly between 0 and 1 (default 0.95).
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
#'     \item{`alpha_ci`}{The lower/upper envelope of grid-supported profile intervals
#'       within `qchisq(conf_level, 1) / 2` log-likelihood units of the maximum.
#'       It includes the refined estimate. Inspect components and truncation flags
#'       before treating the envelope as a bounded, connected interval.}
#'     \item{`alpha_ci_boot`}{The bootstrap percentile interval for `alpha`, or `NULL`
#'       if `n_boot = 0`.}
#'     \item{`boot`}{The vector of bootstrap `alpha` estimates, or `NULL`.}
#'     \item{`logLik`, `AIC`}{The maximized log-likelihood and the AIC of the best model.
#'       `AIC` counts one extra parameter for the estimated `alpha`, so it is comparable
#'       across [ifc()] models and to a null (connectivity-free) model.}
#'     \item{`model`}{The selected engine's fitted model at `alpha`. Treat its coefficients and
#'       standard errors as conditional on `alpha` (see Details).}
#'     \item{`profile`}{A data frame of `alpha` and `logLik` over the grid, used by
#'       `plot()`.}
#'     \item{`boundary`}{`TRUE` if `alpha` sits at a bound of `alpha_range`, a sign the
#'       range should be widened.}
#'     \item{`df`, `nobs`}{The model parameter count plus alpha, and the actual
#'       likelihood sample size (modeled sites for supported unmarked models).
#'       Missing responses do not remove sites from the connectivity network.}
#'     \item{`engine`, `family`}{The fitting engine label and GLM family name;
#'       `family` is `NULL` for unmarked engines.}
#'     \item{`ifc`}{A list with the incidence model number and area transform.}
#'     \item{`alpha_range`, `conf_level`}{The search bounds in distance units and
#'       the requested confidence level.}
#'     \item{`alpha_ci_intervals`}{A data frame of lower/upper endpoints for each
#'       connected set of supported grid points. Separate rows indicate disjoint
#'       support; the reported envelope can span unsupported values.}
#'     \item{`alpha_ci_truncated`}{Named logical lower/upper flags. `TRUE` means
#'       support reaches that search limit; widen the range before interpreting
#'       that endpoint as an uncertainty bound.}
#'     \item{`bootstrap`}{Requested, successful and failed replicate counts.
#'       At least two successful replicates are needed for a percentile interval.}
#'     \item{`diagnostics`}{Rejected candidate count and unique warning/error
#'       reasons. Non-converged or non-finite fits do not enter the objective.}
#'     \item{`comparison`}{Response, modeled rows, likelihood family and weights
#'       used to check that candidate information criteria are comparable.}
#'     \item{`call`}{The original optimization call.}

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
#'   is a grid approximation to likelihood support, with search-limit truncation
#'   and disjoint components reported explicitly. The parametric bootstrap
#'   (`n_boot > 0`) simulates new responses from the fitted model and re-optimizes
#'   `alpha` for each, giving a percentile interval. The bootstrap is parametric on
#'   purpose: connectivity is a whole-network quantity, so resampling sites would break
#'   the network (duplicated sites sit at zero distance). Holding the network fixed and
#'   resampling only the stochastic response preserves that network. GLM/NB
#'   bootstrapping requires a simple response column, such as `y ~ connectivity`.
#'   Transformed or structured response expressions are rejected for bootstrap
#'   rather than transformed a second time. Original missingness and subset
#'   alignment are retained. Failed replicates are counted and reported.
#'
#'   \strong{Caveats.} The GLM treats connectivity as a fixed predictor, so the GLM's own
#'   standard errors and p-values for the `connectivity` coefficient are conditional on
#'   the estimated `alpha`; they do not include uncertainty from estimating alpha. If the response and
#'   `occ_sites` are the same occupancy data, you are predicting occupancy partly from
#'   itself; avoid that.
#'
#'   \strong{Model classes.} The formula method fits a GLM. To use another class, fit it
#'   once with a placeholder `connectivity` term and pass the fitted object: a `glm`, a
#'   negative binomial model from `MASS::glm.nb` (class `negbin`), or supported
#'   single-season `unmarked::occu()` and `unmarked::pcount()` models
#'   (for `unmarked`, put `connectivity` in the `siteCovs`). At each `alpha` the template
#'   is refit from its original call (via [stats::update()] for unmarked), so its
#'   own parameters (the negative binomial theta,
#'   the detection parameters of an occupancy model) are re-estimated each step. Parameter
#'   counts in the information criteria include those parameters plus `alpha`. The
#'   original formula, weights, offsets, subset and controls are retained for GLM/NB
#'   templates. Original data and any objects referenced by the call must remain
#'   available. Unsupported unmarked classes are rejected explicitly.
#'
#'   \strong{Comparing models.} Passing several models (or both area transforms) fits the
#'   best `alpha` for each and ranks them in a model selection table. The table counts the
#'   number of parameters (`K`) correctly: the GLM coefficients, the dispersion parameter
#'   for families that estimate one (such as `gaussian`), and one more for the optimized
#'   `alpha`. Ranking is by AICc, with delta values and normalized weights. If all
#'   candidate AICc values are undefined, selection stops; fit candidates separately
#'   and use a defined criterion with [ifc_modsel()] instead. The same parameter
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
#'                    n_boot = 20)
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
  if(length(formula) != 3L || !("connectivity" %in% all.vars(formula[[3]]))){
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
  .glm_adapter(formula, data, family = family)
}

# Optimize alpha for one (model, area_transform) combination using an adapter.
optimize_engine <- function(adapter, dist_mat, model, area_transform,
                            patch_area, occ_sites, alpha_range,
                            n_grid, n_boot, conf_level) {

  adapter$reset_diagnostics()

  conn_at <- function(alpha) as.numeric(ifc(alpha, dist_mat, model = model,
                     patch_area = patch_area, occ_sites = occ_sites, area_transform = area_transform))
  ll_at <- function(alpha, response = NULL) {
    fit <- adapter$fit(conn_at(alpha), response = response)
    if (is.null(fit)) return(-Inf)
    ll <- adapter$loglik(fit)
    if (length(ll) != 1L || !is.finite(ll)) -Inf else ll
  }
  optim_alpha <- function(response = NULL) {
    grid <- exp(seq(log(alpha_range[1]), log(alpha_range[2]), length.out = n_grid))
    gll <- vapply(grid, function(a) ll_at(a, response), numeric(1))
    if (!any(is.finite(gll))) {
      stop("No candidate model converged with a finite likelihood. ",
           paste(adapter$diagnostics()$reasons, collapse = " "), call. = FALSE)
    }
    k <- which.max(gll)
    opt <- stats::optimize(function(a) {
      ll <- ll_at(a, response)
      if (is.finite(ll)) -ll else 1e100
    }, lower = grid[max(1, k - 1)], upper = grid[min(n_grid, k + 1)])
    # Retain a grid endpoint if refinement finds a worse or invalid candidate.
    candidates <- c(grid[k], opt$minimum)
    scores <- c(gll[k], ll_at(opt$minimum, response))
    best <- which.max(scores)
    list(alpha = candidates[best], logLik = scores[best], grid = grid, gll = gll)
  }
  if (n_boot > 0 && !isTRUE(adapter$simple_response)) {
    stop("Bootstrap requires a simple response column; transformed and structured responses are unsupported.")
  }
  res <- optim_alpha()
  best <- adapter$fit(conn_at(res$alpha))
  if (is.null(best)) stop("The selected candidate failed its convergence check on refitting.")
  thresh <- res$logLik - stats::qchisq(conf_level, 1) / 2
  # Keep disjoint supported intervals visible instead of hiding gaps in one range.
  pr <- data.frame(alpha = c(res$grid, res$alpha), logLik = c(res$gll, res$logLik))
  pr <- pr[order(pr$alpha), ]
  pr <- pr[!duplicated(pr$alpha), ]
  supported <- is.finite(pr$logLik) & pr$logLik >= thresh
  groups <- cumsum(c(TRUE, diff(supported) != 0))
  pieces <- split(pr$alpha[supported], groups[supported])
  intervals <- do.call(rbind, lapply(pieces, function(x) data.frame(lower = min(x), upper = max(x))))
  rownames(intervals) <- NULL
  alpha_ci <- range(pr$alpha[supported])
  truncated <- c(lower = supported[1], upper = supported[nrow(pr)])
  boundary <- res$alpha <= alpha_range[1] * 1.01 || res$alpha >= alpha_range[2] * 0.99
  df <- adapter$npar(best) + 1L
  boot <- NULL
  alpha_ci_boot <- NULL
  bootstrap <- c(requested = n_boot, successful = 0L, failed = 0L)
  if (n_boot > 0) {
    sims <- tryCatch(adapter$simulate(best, n_boot), error = function(e) {
      stop("Cannot simulate an aligned bootstrap response: ", conditionMessage(e), call. = FALSE)
    })
    boot <- vapply(seq_len(n_boot), function(b) {
      tryCatch(optim_alpha(response = sims[[b]])$alpha, error = function(e) NA_real_)
    }, numeric(1))
    bootstrap["successful"] <- sum(is.finite(boot))
    bootstrap["failed"] <- sum(!is.finite(boot))
    if (bootstrap["successful"] >= 2L) {
      a_lo <- (1 - conf_level) / 2
      alpha_ci_boot <- stats::quantile(boot, c(a_lo, 1 - a_lo), na.rm = TRUE)
    } else alpha_ci_boot <- c(NA_real_, NA_real_)
    if (bootstrap["failed"] > 0 || bootstrap["successful"] < 2) warning(
      "Bootstrap: ", bootstrap["successful"], " successful and ", bootstrap["failed"],
      " failed replicates; at least two successful replicates are needed for an interval.", call. = FALSE)
  }
  diagnostics <- adapter$diagnostics()
  if (diagnostics$rejected > 0) warning("Rejected ", diagnostics$rejected,
      " unreliable candidate fits; inspect `$diagnostics`.", call. = FALSE)
  structure(list(alpha = res$alpha, alpha_ci = alpha_ci, alpha_ci_boot = alpha_ci_boot,
                 boot = boot, logLik = res$logLik, df = df, nobs = adapter$nobs(best),
                 AIC = -2 * res$logLik + 2 * df, model = best,
                 profile = data.frame(alpha = res$grid, logLik = res$gll),
                 engine = adapter$label, family = adapter$family,
                 ifc = list(model = model, area_transform = area_transform),
                 alpha_range = alpha_range, conf_level = conf_level, boundary = boundary,
                 alpha_ci_intervals = intervals, alpha_ci_truncated = truncated,
                 bootstrap = bootstrap, diagnostics = diagnostics,
                 comparison = adapter$comparison(best), call = NULL), class = "ifc_optim")
}

# Run the (possibly multi-model) optimization for any adapter, returning a single
# ifc_optim or, for several combinations, an ifc_optim_set.
run_ifc_optim <- function(adapter, dist_mat, model, patch_area, occ_sites,
                          area_transform, alpha_range, n_grid, n_boot,
                          conf_level, call) {

  .check_distances(dist_mat)
  if (adapter$n_sites != nrow(dist_mat)) stop("The model must have the same number of sites as rows in `dist_mat`.")
  if (!is.numeric(model) || !length(model) || any(!is.finite(model)) ||
      any(model != floor(model)) || any(!model %in% 1:4)) stop("`model` must be integer value(s) in 1:4.")
  model <- as.integer(model)
  if (any(model > 1) && is.null(patch_area)) stop("`patch_area` is required for models 2, 3, and 4.")
  .check_scalar(n_grid, "n_grid", lower = 3, integer = TRUE)
  .check_scalar(n_boot, "n_boot", lower = 0, integer = TRUE)
  .check_scalar(conf_level, "conf_level", lower = 0, upper = 1, open_lower = TRUE, open_upper = TRUE)
  area_transform <- match.arg(area_transform, c("none", "log"), several.ok = TRUE)
  dvec <- lower(dist_mat)
  dpos <- dvec[dvec > 0]
  if (!length(dpos)) stop("At least one positive pairwise distance is required to identify a scale.")
  if (is.null(alpha_range)) alpha_range <- c(min(dpos) / 5, max(dvec) * 2)
  if (!is.numeric(alpha_range) || length(alpha_range) != 2L || any(!is.finite(alpha_range)) ||
      alpha_range[1] <= 0 || alpha_range[1] >= alpha_range[2]) {
    stop("`alpha_range` must be a length-2 vector of finite positive, increasing values.")
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
