make_inputs <- function() {
  set.seed(123)
  loc <- data.frame(x = runif(10, 0, 10), y = runif(10, 0, 10))
  list(
    d = dist_mat(loc),
    alpha = 3,
    p_a = runif(10, 5, 50),
    occ = c(0, 0, 1, 1, 0, 1, 0, 0, 1, 1)
  )
}

test_that("ifc returns an ifc object of the right length", {
  x <- make_inputs()
  c1 <- ifc(x$alpha, x$d, model = 1)
  expect_s3_class(c1, "ifc")
  expect_type(as.numeric(c1), "double")
  expect_length(c1, nrow(x$d))
  expect_identical(attr(c1, "model"), 1)
})

test_that("model 1 is the row sums of the connectivity probabilities", {
  x <- make_inputs()
  c_prob <- exp(-(1 / x$alpha) * x$d)
  diag(c_prob) <- 0
  expect_equal(as.numeric(ifc(x$alpha, x$d, model = 1)),
               unname(rowSums(c_prob)))
})

test_that("scaling puts the maximum at 1", {
  x <- make_inputs()
  cs <- ifc(x$alpha, x$d, model = 1, scale = TRUE)
  expect_equal(max(as.numeric(cs)), 1)
  expect_true(attr(cs, "scaled"))
})

test_that("area_transform changes models 2-4 but not model 1", {
  x <- make_inputs()
  m1_none <- as.numeric(ifc(x$alpha, x$d, model = 1))
  m1_log  <- as.numeric(ifc(x$alpha, x$d, model = 1, area_transform = "log"))
  expect_equal(m1_none, m1_log)

  m2_none <- as.numeric(ifc(x$alpha, x$d, model = 2, patch_area = x$p_a))
  m2_log  <- as.numeric(ifc(x$alpha, x$d, model = 2, patch_area = x$p_a,
                            area_transform = "log"))
  expect_false(isTRUE(all.equal(m2_none, m2_log)))
})

test_that("default model 2 uses area directly, not logged", {
  x <- make_inputs()
  c_prob <- exp(-(1 / x$alpha) * x$d)
  diag(c_prob) <- 0
  expected <- rowSums(sweep(c_prob, 2, x$p_a, "*"))
  expect_equal(as.numeric(ifc(x$alpha, x$d, model = 2, patch_area = x$p_a)),
               unname(expected))
})

test_that("occupied sites restrict contributors", {
  x <- make_inputs()
  full <- ifc(x$alpha, x$d, model = 1)
  occ  <- ifc(x$alpha, x$d, model = 1, occ_sites = x$occ)
  expect_false(isTRUE(all.equal(as.numeric(full), as.numeric(occ))))
})

test_that("ifc validates its inputs", {
  x <- make_inputs()
  expect_error(ifc(x$alpha, as.vector(x$d), model = 1), "square")
  expect_error(ifc(x$alpha, x$d, model = 5), "1-4")
  expect_error(ifc(x$alpha, x$d, model = 2), "patch_area")
  expect_error(ifc(x$alpha, x$d, model = 1, occ_sites = c(1, 0)), "length")
  expect_error(ifc(x$alpha, x$d, model = 2, patch_area = c(1, 2)), "length")
})

test_that("ifc methods run without error", {
  x <- make_inputs()
  c4 <- ifc(x$alpha, x$d, model = 4, patch_area = x$p_a)
  expect_output(print(c4), "Incidence function connectivity")
  expect_output(summary(c4), "summary")
  expect_invisible(plot(c4))
})
