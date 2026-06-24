sim_set <- function(n = 120, true_alpha = 4, seed = 1) {
  set.seed(seed)
  loc  <- data.frame(x = runif(n, 0, 20), y = runif(n, 0, 20))
  d    <- dist_mat(loc)
  conn <- as.numeric(ifc(true_alpha, d, model = 1, scale = TRUE))
  dat  <- data.frame(abundance = 2 + 8 * conn + rnorm(n, 0, 0.5))
  list(d = d, dat = dat, p_a = runif(n, 5, 50))
}

test_that("logLik counts the alpha parameter", {
  s <- sim_set()
  fit <- ifc_optim(abundance ~ connectivity, data = s$dat, dist_mat = s$d,
                   model = 1, family = gaussian)
  expect_equal(attr(logLik(fit), "df"), attr(logLik(fit$model), "df") + 1L)
  expect_equal(nobs(fit), nobs(fit$model))
  expect_equal(AIC(fit), AIC(fit$model) + 2)        # +1 parameter -> +2 AIC
  expect_equal(fit$AIC, AIC(fit))                   # stored value agrees
})

test_that("AICc default works and exceeds AIC", {
  s <- sim_set()
  fit <- ifc_optim(abundance ~ connectivity, data = s$dat, dist_mat = s$d,
                   model = 1, family = gaussian)
  k <- attr(logLik(fit), "df"); n <- nobs(fit)
  expect_equal(AICc(fit), AIC(fit) + (2 * k * (k + 1)) / (n - k - 1))
  expect_gt(AICc(fit), AIC(fit))
  ## default method also works on a plain glm
  g <- glm(abundance ~ 1, data = s$dat, family = gaussian)
  expect_true(is.finite(AICc(g)))
})

test_that("a vector of models returns a ranked ifc_optim_set", {
  s <- sim_set()
  set <- ifc_optim(abundance ~ connectivity, data = s$dat, dist_mat = s$d,
                   model = 1:4, patch_area = s$p_a, family = gaussian)
  expect_s3_class(set, "ifc_optim_set")
  expect_length(set$fits, 4)                        # model 1 + models 2,3,4 (none)
  expect_s3_class(set$table, "ifc_modsel")
  expect_s3_class(set$best, "ifc_optim")
  ## ranked: best row has delta 0 and the largest weight
  expect_equal(set$table$delta[1], 0)
  expect_equal(which.max(set$table$weight), 1L)
  expect_equal(sum(set$table$weight), 1, tolerance = 1e-8)
  expect_identical(set$best, set$fits[[set$table$model[1]]])
})

test_that("both area transforms expand the comparison", {
  s <- sim_set()
  set <- ifc_optim(abundance ~ connectivity, data = s$dat, dist_mat = s$d,
                   model = 1:4, patch_area = s$p_a,
                   area_transform = c("none", "log"), family = gaussian)
  ## model 1 once + models 2,3,4 x 2 transforms = 7
  expect_length(set$fits, 7)
})

test_that("ifc_modsel ranks loose fits and counts K with alpha", {
  s <- sim_set()
  m1 <- ifc_optim(abundance ~ connectivity, data = s$dat, dist_mat = s$d,
                  model = 1, family = gaussian)
  m2 <- ifc_optim(abundance ~ connectivity, data = s$dat, dist_mat = s$d,
                  model = 2, patch_area = s$p_a, family = gaussian)
  tab <- ifc_modsel(m1, m2)
  expect_s3_class(tab, "ifc_modsel")
  expect_equal(nrow(tab), 2)
  expect_true(all(tab$K == 4))                      # 2 coefs + dispersion + alpha
  expect_equal(min(tab$delta), 0)
})

test_that("ifc_optim errors when models 2-4 lack patch_area", {
  s <- sim_set()
  expect_error(
    ifc_optim(abundance ~ connectivity, data = s$dat, dist_mat = s$d, model = 1:4,
              family = gaussian),
    "patch_area")
})

test_that("ifc_optim_set methods run without error", {
  s <- sim_set()
  set <- ifc_optim(abundance ~ connectivity, data = s$dat, dist_mat = s$d,
                   model = 1:4, patch_area = s$p_a, family = gaussian)
  expect_output(print(set), "model selection")
  expect_output(summary(set), "Best model")
  expect_invisible(plot(set))
})
