# Tests for the fitted-model template engines (glm, glm.nb, unmarked).

make_sites <- function(n = 100, true_alpha = 4, seed = 1) {
  set.seed(seed)
  loc  <- data.frame(x = runif(n, 0, 20), y = runif(n, 0, 20))
  d    <- dist_mat(loc)
  conn <- as.numeric(ifc(true_alpha, d, model = 1, scale = TRUE))
  list(d = d, conn = conn, n = n)
}

test_that("a glm template gives the same result as the formula method", {
  s <- make_sites()
  dat <- data.frame(abundance = 2 + 8 * s$conn + rnorm(s$n, 0, 0.5))
  f_formula <- ifc_optim(abundance ~ connectivity, data = dat, dist_mat = s$d,
                         model = 1, family = gaussian)

  dat0 <- dat; dat0$connectivity <- 0
  g0   <- glm(abundance ~ connectivity, data = dat0, family = gaussian)
  f_obj <- ifc_optim(g0, dist_mat = s$d, model = 1)

  expect_s3_class(f_obj, "ifc_optim")
  expect_equal(f_obj$alpha, f_formula$alpha, tolerance = 1e-4)
  expect_equal(f_obj$df, f_formula$df)
})

test_that("ifc_optim errors on a template without a connectivity term", {
  s <- make_sites()
  dat <- data.frame(y = rnorm(s$n), x = rnorm(s$n))
  g <- glm(y ~ x, data = dat, family = gaussian)
  expect_error(ifc_optim(g, dist_mat = s$d), "connectivity")
})

test_that("ifc_optim errors for an unsupported class", {
  s <- make_sites()
  expect_error(ifc_optim(1:5, dist_mat = s$d), "No `ifc_optim\\(\\)` method")
})

test_that("glm.nb template counts theta and alpha", {
  skip_if_not_installed("MASS")
  s <- make_sites()
  mu <- exp(0.5 + 3 * s$conn)
  dn <- data.frame(count = rnbinom(s$n, mu = mu, size = 3), connectivity = 0)
  nb0 <- MASS::glm.nb(count ~ connectivity, data = dn)

  fit <- ifc_optim(nb0, dist_mat = s$d, model = 1)
  expect_s3_class(fit, "ifc_optim")
  expect_identical(fit$engine, "negative binomial")
  ## 1 intercept + 1 connectivity + theta + alpha = 4
  expect_equal(fit$df, 4L)
  expect_equal(AIC(fit), -2 * fit$logLik + 2 * fit$df)
  expect_true(is.finite(fit$alpha) && fit$alpha > 0)
})

test_that("unmarked occu template is supported", {
  skip_if_not_installed("unmarked")
  s <- make_sites()
  set.seed(7)
  z <- rbinom(s$n, 1, plogis(-0.5 + 4 * s$conn))
  y <- matrix(NA, s$n, 3)
  for (j in 1:3) y[, j] <- z * rbinom(s$n, 1, 0.6)
  umf <- unmarked::unmarkedFrameOccu(
    y = y, siteCovs = data.frame(connectivity = numeric(s$n)))
  ## The placeholder connectivity is constant, so the template fit is degenerate
  ## (singular Hessian); ifc_optim() refits with real connectivity. Silence it here.
  oc0 <- suppressWarnings(unmarked::occu(~1 ~ connectivity, data = umf))

  fit <- ifc_optim(oc0, dist_mat = s$d, model = 1)
  expect_s3_class(fit, "ifc_optim")
  expect_identical(fit$engine, "unmarkedFitOccu")
  expect_equal(fit$nobs, s$n)
  ## 2 occupancy params + 1 detection param + alpha = 4
  expect_equal(fit$df, 4L)
  expect_true(is.finite(as.numeric(logLik(fit))))
  expect_equal(AIC(fit), -2 * fit$logLik + 2 * fit$df)
})

test_that("multi-model selection works through a template", {
  skip_if_not_installed("MASS")
  s <- make_sites()
  mu <- exp(0.5 + 3 * s$conn)
  dn <- data.frame(count = rnbinom(s$n, mu = mu, size = 3), connectivity = 0)
  nb0 <- MASS::glm.nb(count ~ connectivity, data = dn)

  set <- ifc_optim(nb0, dist_mat = s$d, model = 1:4, patch_area = runif(s$n, 5, 50))
  expect_s3_class(set, "ifc_optim_set")
  expect_length(set$fits, 4)
  expect_s3_class(set$table, "ifc_modsel")
})
