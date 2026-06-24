test_that("lower extracts the lower triangle", {
  m <- matrix(c(0, 1, 2,
                1, 0, 3,
                2, 3, 0), 3, 3, byrow = TRUE)
  expect_equal(sort(lower(m)), c(1, 2, 3))
  expect_length(lower(m), choose(3, 2))
})

test_that("lower rejects a non-square matrix (regression for the [1] vs [2] bug)", {
  m <- matrix(1:6, nrow = 2, ncol = 3)
  expect_error(lower(m), "square")
})
