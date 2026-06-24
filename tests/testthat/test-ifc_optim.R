## A continuous response identifies the dispersal scale well, so these tests are
## stable. (Binary occupancy carries little information about alpha; that flat-profile
## case is intentional behavior, documented in ?ifc_optim, not tested for recovery.)
sim_abundance <- function(n = 120, true_alpha = 4, seed = 1) {
  set.seed(seed)
  loc  <- data.frame(x = runif(n, 0, 20), y = runif(n, 0, 20))
  d    <- dist_mat(loc)
  conn <- as.numeric(ifc(true_alpha, d, model = 1, scale = TRUE))
  dat  <- data.frame(abundance = 2 + 8 * conn + rnorm(n, 0, 0.5))
  list(d = d, dat = dat, true_alpha = true_alpha)
}

test_that("ifc_optim returns a well-formed ifc_optim object", {
  s <- sim_abundance()
  fit <- ifc_optim(abundance ~ connectivity, data = s$dat, dist_mat = s$d,
                   model = 1, family = gaussian)
  expect_s3_class(fit, "ifc_optim")
  expect_true(is.finite(fit$alpha))
  expect_gt(fit$alpha, 0)
  expect_s3_class(fit$model, "glm")
  expect_length(fit$alpha_ci, 2)
})

test_that("the optimized alpha beats the range endpoints in log-likelihood", {
  s <- sim_abundance()
  fit <- ifc_optim(abundance ~ connectivity, data = s$dat, dist_mat = s$d,
                   model = 1, family = gaussian)
  endpoints <- fit$profile$logLik[c(1, nrow(fit$profile))]
  expect_true(fit$logLik >= max(endpoints, na.rm = TRUE))
})

test_that("alpha is recovered with an informative response", {
  s <- sim_abundance(n = 120, true_alpha = 4)
  fit <- ifc_optim(abundance ~ connectivity, data = s$dat, dist_mat = s$d,
                   model = 1, family = gaussian)
  ## generous bound: within a factor of 2 of the truth
  expect_gt(fit$alpha, s$true_alpha / 2)
  expect_lt(fit$alpha, s$true_alpha * 2)
})

test_that("AIC counts the alpha parameter", {
  s <- sim_abundance()
  fit <- ifc_optim(abundance ~ connectivity, data = s$dat, dist_mat = s$d,
                   model = 1, family = gaussian)
  expect_equal(fit$AIC, AIC(fit$model) + 2)
})

test_that("bootstrap produces a CI and replicate vector", {
  s <- sim_abundance()
  fit <- ifc_optim(abundance ~ connectivity, data = s$dat, dist_mat = s$d,
                   model = 1, family = gaussian, n_boot = 25)
  expect_length(fit$boot, 25)
  expect_length(fit$alpha_ci_boot, 2)
  expect_true(fit$alpha_ci_boot[1] <= fit$alpha_ci_boot[2])
})

test_that("ifc_optim validates its inputs", {
  s <- sim_abundance()
  expect_error(ifc_optim(abundance ~ habitat, data = s$dat, dist_mat = s$d),
               "connectivity")
  expect_error(ifc_optim(abundance ~ connectivity, data = s$dat[1:5, , drop = FALSE],
                         dist_mat = s$d), "same number of rows")
  expect_error(ifc_optim(abundance ~ connectivity, data = s$dat, dist_mat = s$d,
                         alpha_range = c(5, 1)), "increasing")
})

test_that("ifc_optim methods run without error", {
  s <- sim_abundance()
  fit <- ifc_optim(abundance ~ connectivity, data = s$dat, dist_mat = s$d,
                   model = 1, family = gaussian)
  expect_output(print(fit), "alpha")
  expect_output(summary(fit), "GLM")
  expect_invisible(plot(fit))
})
