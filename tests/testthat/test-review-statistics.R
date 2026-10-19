review_data <- function(n = 60) {
  set.seed(42)
  d <- dist_mat(cbind(runif(n, 0, 20), runif(n, 0, 20)))
  conn <- as.numeric(ifc(3, d, model = 1))
  dat <- data.frame(y = 2 + .2 * conn + rnorm(n, sd = .3), connectivity = conn,
                    w = seq(.2, 3, length.out = n), off = rnorm(n, sd = .2), keep = seq_len(n) > 4)
  list(d = d, dat = dat)
}

test_that("scientific input domains are enforced and zero connectivity scales to zero", {
  d <- as.matrix(dist(cbind(1:3, 3:1)))
  for (a in list(0, -1, Inf, NA_real_, c(1, 2))) expect_error(ifc(a, d, 1), "alpha")
  expect_error(ifc(2, d, 1.9), "integer")
  expect_error(ifc(2, d, 1, occ_sites = c(-1, 1, 1)), "0 or 1")
  expect_error(ifc(2, d, 2, patch_area = c(1, -1, 2)), "positive")
  expect_equal(as.numeric(ifc(2, d, 1, occ_sites = rep(0, 3), scale = TRUE)), rep(0, 3))
  for (bad in list(data.frame(x = letters[1:3], y = 1:3), cbind(c(1, NA, 3), 1:3))) {
    expect_error(dist_mat(bad), "finite numeric")
  }
  expect_error(lower(1:3), "square")
  bad <- d; bad[1, 2] <- -1
  expect_error(ifc(2, bad, 1), "nonnegative")
  expect_error(ifc(2, diag(3), 1), "zero diagonal")
})

test_that("template weights, offsets, subsets, and controls survive refitting", {
  x <- review_data()
  g <- glm(y ~ connectivity, data = x$dat, weights = w, offset = off,
            subset = keep, family = gaussian, control = glm.control(epsilon = 1e-9, maxit = 40))
  f <- ifc_optim(g, dist_mat = x$d, alpha_range = c(.5, 10), n_grid = 12)
  ref <- x$dat
  ref$connectivity <- as.numeric(ifc(f$alpha, x$d, 1))
  direct <- glm(y ~ connectivity, data = ref, weights = w, offset = off,
                 subset = keep, family = gaussian, control = glm.control(epsilon = 1e-9, maxit = 40))
  expect_equal(f$model$prior.weights, direct$prior.weights)
  expect_equal(f$model$offset, direct$offset)
  expect_equal(coef(f$model), coef(direct))
  expect_equal(f$model$control, direct$control)
  expect_equal(f$nobs, 56)
})

test_that("missing responses retain their network rows during bootstrap", {
  x <- review_data()
  x$dat$y[1:5] <- NA
  f <- ifc_optim(y ~ connectivity, data = x$dat, dist_mat = x$d, family = gaussian,
                  alpha_range = c(.5, 10), n_grid = 12, n_boot = 4)
  expect_equal(nobs(f), nobs(f$model))
  expect_equal(nobs(f), 55)
  expect_equal(BIC(f), -2 * f$logLik + f$df * log(55))
  expect_equal(unname(f$bootstrap), c(4, 4, 0))
  expect_true(all(is.finite(f$boot)))
  a <- adapter_formula(y ~ connectivity, x$dat, gaussian())
  sims <- a$simulate(f$model, 2)
  expect_equal(length(sims[[1]]), 60)
  expect_true(all(is.na(sims[[1]][1:5])))
})

test_that("unsupported bootstrap response expressions and invalid controls fail early", {
  x <- review_data()
  x$dat$y <- exp(x$dat$y)
  expect_error(ifc_optim(log(y) ~ connectivity, data = x$dat, dist_mat = x$d,
                          family = gaussian, n_boot = 2), "simple response column")
  for (args in list(list(model = 1.9), list(n_grid = 1), list(conf_level = 2), list(n_boot = -1))) {
    expect_error(do.call(ifc_optim, c(list(object = y ~ connectivity, data = x$dat,
                                         dist_mat = x$d, family = gaussian), args)), "must be")
  }
})

test_that("separated nonconvergent candidates cannot become selected models", {
  set.seed(1)
  d <- dist_mat(cbind(runif(60), runif(60)))
  conn <- as.numeric(ifc(.15, d, 1))
  dat <- data.frame(y = as.integer(conn > median(conn)), connectivity = conn)
  # Deliberately insufficient iterations give a reliable convergence regression.
  g <- suppressWarnings(glm(y ~ connectivity, data = dat, family = binomial,
                             control = glm.control(maxit = 1)))
  expect_error(ifc_optim(g, dist_mat = d, alpha_range = c(.01, 1), n_grid = 12),
                "No candidate model converged")
})

test_that("profile interval components and range truncation are explicit", {
  x <- review_data()
  f <- ifc_optim(y ~ connectivity, data = x$dat, dist_mat = x$d, family = gaussian,
                  alpha_range = c(2.9, 3.1), n_grid = 8)
  expect_equal(names(f$alpha_ci_intervals), c("lower", "upper"))
  expect_named(f$alpha_ci_truncated, c("lower", "upper"))
  expect_true(any(f$alpha_ci_truncated))
  expect_true(f$model$converged)
  expect_equal(f$diagnostics$rejected, 0)
})

test_that("model selection refuses incomparable responses and undefined rankings", {
  x <- review_data()
  f <- ifc_optim(y ~ connectivity, data = x$dat, dist_mat = x$d, family = gaussian, n_grid = 12)
  y <- x$dat; y$y <- 2 * y$y
  g <- ifc_optim(y ~ connectivity, data = y, dist_mat = x$d, family = gaussian, n_grid = 12)
  expect_error(ifc_modsel(f, g), "same response")
  old <- f; old$comparison <- NULL
  expect_error(ifc_modsel(old), "Refit older")
  old$alpha_ci_intervals <- old$alpha_ci_truncated <- old$diagnostics <- old$bootstrap <- NULL
  old$boot <- c(2, 3, NA_real_); old$alpha_ci_boot <- c(2, 3)
  expect_output(print(old), "2 successful; 1 failed")
  tiny <- f; tiny$nobs <- 3L
  expect_error(ifc_modsel(tiny), "No finite AICc")
  expect_true(is.finite(ifc_modsel(tiny, rank = "AIC")$AIC))
})

test_that("unmarked validates the state formula and preserves missing observations", {
  skip_if_not_installed("unmarked")
  x <- review_data(50)
  set.seed(3)
  y <- matrix(rbinom(150, 1, .4), 50, 3)
  y[1, 1] <- NA
  umf <- unmarked::unmarkedFrameOccu(y, siteCovs = data.frame(connectivity = as.numeric(ifc(3, x$d, 1))))
  unused <- unmarked::occu(~1 ~1, data = umf)
  expect_error(ifc_optim(unused, dist_mat = x$d), "state formula")
  used <- unmarked::occu(~1 ~connectivity, data = umf)
  a <- adapter_unmarked(used, umf)
  simulated <- y; simulated[1, 1] <- 1
  z <- a$fit(as.numeric(ifc(3, x$d, 1)), response = simulated)
  expect_true(is.na(unmarked::getY(unmarked::getData(z))[1, 1]))
})

test_that("negative binomial templates retain offsets and weights", {
  skip_if_not_installed("MASS")
  x <- review_data(100)
  set.seed(7)
  x$dat$y <- rnbinom(100, mu = exp(.3 + .15 * x$dat$connectivity + x$dat$off), size = 3)
  g <- MASS::glm.nb(y ~ connectivity + offset(off), data = x$dat, weights = w)
  f <- ifc_optim(g, dist_mat = x$d, n_grid = 12, alpha_range = c(1, 8))
  expect_equal(f$model$prior.weights, g$prior.weights)
  expect_equal(f$model$offset, g$offset)
  expect_equal(f$df, 4)
})

test_that("supported unmarked engines bootstrap aligned response matrices", {
  skip_if_not_installed("unmarked")
  x <- review_data(60)
  conn <- as.numeric(ifc(3, x$d, 1, scale = TRUE))
  set.seed(12)
  n <- rpois(60, exp(.3 + conn))
  counts <- sapply(1:3, function(j) rbinom(60, n, .7))
  detected <- (counts > 0) + 0
  counts[1, 1] <- detected[1, 1] <- NA
  p <- unmarked::unmarkedFramePCount(counts, siteCovs = data.frame(connectivity = conn))
  o <- unmarked::unmarkedFrameOccu(detected, siteCovs = data.frame(connectivity = conn))
  templates <- list(unmarked::pcount(~1 ~connectivity, data = p, K = 25),
                     unmarked::occu(~1 ~connectivity, data = o))
  for (template in templates) {
    messages <- character()
    f <- withCallingHandlers(ifc_optim(template, dist_mat = x$d, alpha_range = c(1, 8),
                                        n_grid = 8, n_boot = 2), warning = function(w) {
      messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
    })
    expect_equal(unname(f$bootstrap), c(2, 2, 0))
    expect_equal(f$nobs, 60)
    expect_true(all(is.finite(f$boot)))
    expect_equal(f$model@opt$convergence, 0)
    if (f$diagnostics$rejected > 0) expect_true(any(grepl("unreliable candidate", messages)))
    else expect_length(messages, 0)
  }
})
