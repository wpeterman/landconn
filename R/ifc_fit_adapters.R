# Catch candidate-specific failures without concealing them from the final result.
# No fit with a failed convergence check contributes to the objective or interval.
.guard_fits <- function(fit, valid) {
  state <- new.env(parent = emptyenv())
  state$rejected <- 0L
  state$reasons <- character()
  guarded <- function(conn, response = NULL) {
    warnings <- character()
    reason <- NULL
    ans <- tryCatch(withCallingHandlers(fit(conn, response), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) { reason <<- conditionMessage(e); NULL })
    separation <- any(grepl("did not converge|numerically 0 or 1|iteration limit|alternation limit", warnings))
    if (is.null(ans) || separation || !isTRUE(valid(ans))) {
      state$rejected <- state$rejected + 1L
      state$reasons <- unique(c(state$reasons, reason, warnings, "Failed convergence or finite-likelihood check."))
      return(NULL)
    }
    # Other warnings remain visible in the diagnostics even for usable fits.
    state$reasons <- unique(c(state$reasons, warnings))
    ans
  }
  list(fit = guarded, diagnostics = function() list(rejected = state$rejected, reasons = state$reasons),
       reset_diagnostics = function() { state$rejected <- 0L; state$reasons <- character(); invisible(NULL) })
}

.glm_adapter <- function(formula, data, family = NULL, template = NULL, negative_binomial = FALSE) {
  simple <- is.symbol(formula[[2]])
  response_name <- if (simple) as.character(formula[[2]]) else NULL
  original_response <- if (simple) data[[response_name]] else NULL
  if (negative_binomial && !requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' is required.")
  if (!is.null(template)) {
    call <- stats::getCall(template)
    call[[1]] <- if (negative_binomial) quote(MASS::glm.nb) else quote(stats::glm)
    # Preserve subset, weights, offset, control, link, and other fitted-call settings.
    # Only the data and formula are replaced; expressions evaluate in their original
    # formula environment, with the full site-network data available.
    call$data <- quote(.landconn_data)
    call$formula <- formula
  }
  fit <- function(conn, response = NULL) {
    dat <- data
    dat$connectivity <- conn
    if (!is.null(response)) {
      if (!simple) stop("Bootstrap requires a simple response column; transformed and structured responses are unsupported.")
      if (length(response) != nrow(data)) stop("Simulated response is not aligned with the full site network.")
      dat[[response_name]] <- response
    }
    if (is.null(template)) {
      if (negative_binomial) MASS::glm.nb(formula, data = dat) else stats::glm(formula, data = dat, family = family)
    } else {
      env <- new.env(parent = environment(formula))
      env$.landconn_data <- dat
      eval(call, env)
    }
  }
  valid <- function(fitted) {
    isTRUE(fitted$converged) && !isTRUE(fitted$boundary) && is.finite(as.numeric(stats::logLik(fitted))) &&
      (!negative_binomial || (is.null(fitted$th.warn) && is.finite(fitted$theta) && fitted$theta > 0))
  }
  guarded <- .guard_fits(fit, valid)
  label <- if (negative_binomial) "negative binomial" else paste0("glm (", family$family, ")")
  c(list(n_sites = nrow(data), label = label,
         family = if (negative_binomial) "Negative Binomial" else family$family,
         simple_response = simple,
         loglik = function(fitted) as.numeric(stats::logLik(fitted)),
         npar = function(fitted) as.integer(attr(stats::logLik(fitted), "df")),
         nobs = function(fitted) stats::nobs(fitted),
         comparison = function(fitted) list(engine = if (negative_binomial) "negative binomial" else family$family,
              response = unname(stats::model.response(stats::model.frame(fitted))),
              rows = rownames(stats::model.frame(fitted)), weights = unname(fitted$prior.weights)),
         simulate = function(fitted, nsim) {
           if (!simple) stop("Bootstrap requires a simple response column.")
           sims <- stats::simulate(fitted, nsim = nsim)
           rows <- match(rownames(stats::model.frame(fitted)), rownames(data))
           if (anyNA(rows)) stop("Cannot align fitted rows with the original site network.")
           lapply(seq_len(nsim), function(i) {
             if (!is.null(dim(sims[[i]]))) stop("Bootstrap for structured responses is unsupported.")
             ans <- original_response
             ans[rows] <- sims[[i]]
             ans
           })
         }), guarded)
}
