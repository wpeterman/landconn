# ifc_optim engines: fitted-model template methods -----------------------------

# These methods let ifc_optim() refit an already fitted model across alpha. The user
# fits the model once with a placeholder `connectivity` term; each method builds an
# adapter (see adapter_formula() in ifc_optim.R for the interface) and hands it to
# run_ifc_optim(). See ?ifc_optim.

# Recover the data frame a fitted model was built from. Prefer the original data still
# in scope; fall back to the model frame.
.template_data <- function(object) {
  d <- tryCatch(eval(stats::getCall(object)$data, environment(stats::formula(object))),
                error = function(e) NULL)
  if(is.null(d) || !is.data.frame(d)) d <- stats::model.frame(object)
  d
}

#' @rdname ifc_optim
#' @exportS3Method ifc_optim glm
#' @importFrom stats family formula getCall model.frame update
ifc_optim.glm <- function(object, dist_mat, model = 1, patch_area = NULL,
                          occ_sites = NULL, area_transform = "none",
                          alpha_range = NULL, n_grid = 50, n_boot = 0,
                          conf_level = 0.95, ...) {
  formula <- stats::formula(object)
  if(length(formula) != 3L || !("connectivity" %in% all.vars(formula[[3]]))){
    stop("The fitted model must include a `connectivity` term. Fit it once with a placeholder `connectivity` covariate, then pass it here.")
  }
  data    <- .template_data(object)
  adapter <- .glm_adapter(formula, data, family = stats::family(object), template = object)
  run_ifc_optim(adapter, dist_mat, model = model, patch_area = patch_area,
                occ_sites = occ_sites, area_transform = area_transform,
                alpha_range = alpha_range, n_grid = n_grid, n_boot = n_boot,
                conf_level = conf_level, call = match.call())
}

#' @rdname ifc_optim
#' @exportS3Method ifc_optim negbin
ifc_optim.negbin <- function(object, dist_mat, model = 1, patch_area = NULL,
                             occ_sites = NULL, area_transform = "none",
                             alpha_range = NULL, n_grid = 50, n_boot = 0,
                             conf_level = 0.95, ...) {
  if(!requireNamespace("MASS", quietly = TRUE)){
    stop("Package 'MASS' is required to refit a negative binomial model.")
  }
  formula <- stats::formula(object)
  if(length(formula) != 3L || !("connectivity" %in% all.vars(formula[[3]]))){
    stop("The fitted model must include a `connectivity` term. Fit it once with a placeholder `connectivity` covariate, then pass it here.")
  }
  data    <- .template_data(object)
  adapter <- .glm_adapter(formula, data, template = object, negative_binomial = TRUE)
  run_ifc_optim(adapter, dist_mat, model = model, patch_area = patch_area,
                occ_sites = occ_sites, area_transform = area_transform,
                alpha_range = alpha_range, n_grid = n_grid, n_boot = n_boot,
                conf_level = conf_level, call = match.call())
}

# Adapter that refits via MASS::glm.nb. The negative binomial theta is re-estimated at
# each alpha; logLik() already counts it in the degrees of freedom.
adapter_negbin <- function(formula, data) {
  .glm_adapter(formula, data, negative_binomial = TRUE)
}

#' @rdname ifc_optim
#' @exportS3Method ifc_optim default
ifc_optim.default <- function(object, ...) {
  ## unmarked fits are S4, so they reach the default method; route them by S4 class.
  if(methods::is(object, "unmarkedFit")){
    return(ifc_optim_unmarked(object, ...))
  }
  stop("No `ifc_optim()` method for an object of class ",
       paste(sQuote(class(object)), collapse = "/"),
       ". Supply a formula, or a fitted `glm`, `glm.nb`, or `unmarked` model.")
}

# Handler for any unmarkedFit subclass (occu, pcount, ...). Not an S3 method because
# unmarked fits are S4; ifc_optim.default routes here via methods::is().
ifc_optim_unmarked <- function(object, dist_mat, model = 1, patch_area = NULL,
                               occ_sites = NULL, area_transform = "none",
                               alpha_range = NULL, n_grid = 50, n_boot = 0,
                               conf_level = 0.95, ...) {
  if(!requireNamespace("unmarked", quietly = TRUE)){
    stop("Package 'unmarked' is required for unmarked models.")
  }
  if (!class(object)[1] %in% c("unmarkedFitOccu", "unmarkedFitPCount")) {
    stop("Supported unmarked templates are single-season occu() and pcount() fits.")
  }
  formula <- methods::slot(object, "formula")
  if (length(formula) != 3L || !"connectivity" %in% all.vars(formula[[3]]) ||
      "connectivity" %in% all.vars(formula[[2]])) {
    stop("The unmarked state formula must use `connectivity`; the detection formula must not use it.")
  }
  umf <- unmarked::getData(object)
  sc  <- unmarked::siteCovs(umf)
  if(is.null(sc) || !("connectivity" %in% names(sc))){
    stop("The unmarked model's `siteCovs` must include a `connectivity` column. Fit it once with a placeholder, then pass it here.")
  }
  adapter <- adapter_unmarked(object, umf)
  run_ifc_optim(adapter, dist_mat, model = model, patch_area = patch_area,
                occ_sites = occ_sites, area_transform = area_transform,
                alpha_range = alpha_range, n_grid = n_grid, n_boot = n_boot,
                conf_level = conf_level, call = match.call())
}

# Adapter that refits an unmarked model via update(). Connectivity is a site covariate;
# the bootstrap replaces the response matrix in the unmarkedFrame. Log-likelihood and the
# parameter count are read from S4 slots: `unmarked` defines logLik() and coef() as S4
# methods, which a namespace-qualified `stats::logLik()` would bypass.
adapter_unmarked <- function(object, umf) {
  original_y <- unmarked::getY(umf)
  fit <- function(conn, response = NULL) {
    u <- umf
    sc <- unmarked::siteCovs(u)
    sc$connectivity <- conn
    unmarked::siteCovs(u) <- sc
    if (!is.null(response)) {
      if (!identical(dim(response), dim(original_y))) stop("Simulated unmarked response dimensions changed.")
      response[is.na(original_y)] <- NA
      methods::slot(u, "y") <- response
    }
    stats::update(object, data = u)
  }
  check <- function(fitted) {
    opt <- methods::slot(fitted, "opt")
    isTRUE(opt$convergence == 0) && is.finite(methods::slot(fitted, "negLogLike")) &&
      all(is.finite(opt$par))
  }
  guarded <- .guard_fits(fit, check)
  c(list(n_sites = unmarked::numSites(umf), label = class(object)[1], family = NULL,
         simple_response = TRUE,
         loglik = function(fitted) -as.numeric(methods::slot(fitted, "negLogLike")),
         npar = function(fitted) as.integer(round(methods::slot(fitted, "AIC") / 2 - methods::slot(fitted, "negLogLike"))),
         nobs = function(fitted) unmarked::numSites(unmarked::getData(fitted)) - length(methods::slot(fitted, "sitesRemoved")),
         comparison = function(fitted) list(engine = class(fitted)[1], response = original_y,
              rows = setdiff(seq_len(nrow(original_y)), methods::slot(fitted, "sitesRemoved"))),
         simulate = function(fitted, nsim) unmarked::simulate(fitted, nsim = nsim)), guarded)
}
