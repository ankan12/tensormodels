test_that("basic example", {
  a <- array(1:3, dim = c(3, 1, 1))
  b <- matrix(4:9, nrow = 2, ncol = 3)

  expect_equal(n_mode_prod(a, b, 1),
               array(data = c(40, 46), dim = c(2, 1, 1)))
})
