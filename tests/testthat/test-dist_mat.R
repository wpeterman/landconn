test_that("dist_mat returns a square land_dist matrix", {
  set.seed(1)
  loc <- data.frame(x = runif(8), y = runif(8))
  d <- dist_mat(loc)

  expect_s3_class(d, "land_dist")
  expect_true(inherits(d, "matrix"))
  expect_equal(dim(d), c(8, 8))
  expect_equal(unname(diag(d)), rep(0, 8))
  expect_equal(unclass(d), t(unclass(d)))   # symmetric
})

test_that("dist_mat rejects non two-column input", {
  expect_error(dist_mat(data.frame(x = 1:3)), "two columns")
  expect_error(dist_mat(matrix(1:9, 3, 3)), "two columns")
})

test_that("a land_dist matrix can be passed straight to ifc and lower", {
  set.seed(2)
  loc <- data.frame(x = runif(6), y = runif(6))
  d <- dist_mat(loc)
  expect_s3_class(ifc(1, d, model = 1), "ifc")
  expect_length(lower(d), choose(6, 2))
})

test_that("land_dist methods run without error", {
  set.seed(3)
  loc <- data.frame(x = runif(12), y = runif(12))
  d <- dist_mat(loc)
  expect_output(print(d), "distance matrix")
  expect_output(summary(d), "pairwise distances")
  expect_invisible(plot(d))
})
