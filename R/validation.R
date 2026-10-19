# Shared input contracts keep direct connectivity and fitted workflows aligned.
.check_distances <- function(x) {
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) != ncol(x) || !nrow(x)) {
    stop("`dist_mat` must be a nonempty square numeric distance matrix.", call. = FALSE)
  }
  if (any(!is.finite(x)) || any(x < 0) || any(abs(diag(x)) > 1e-8) ||
      !isTRUE(all.equal(unname(x), unname(t(x)), check.attributes = FALSE, tolerance = 1e-8))) {
    stop("`dist_mat` must be finite, nonnegative, symmetric, and have a zero diagonal.", call. = FALSE)
  }
  invisible(NULL)
}

.check_scalar <- function(x, name, lower = -Inf, upper = Inf, integer = FALSE,
                           open_lower = FALSE, open_upper = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x > upper || (open_lower && x == lower) || (open_upper && x == upper) ||
      (integer && x != floor(x))) {
    stop("`", name, "` must be a finite ", if (integer) "integer" else "number",
         " in the documented range.", call. = FALSE)
  }
  invisible(NULL)
}
