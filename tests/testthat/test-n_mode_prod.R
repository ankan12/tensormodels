test_that("basic example", {
  a <- array(1:3, dim = c(3, 1, 1))
  b <- matrix(4:9, nrow = 2, ncol = 3)

  expect_equal(n_prod(a, 1, b),
               array(data = c(40, 46), dim = c(2, 1, 1)))
})

test_that("formula infix n-mode product", {
  a <- array(1:3, dim = c(3, 1, 1))
  b <- matrix(4:9, nrow = 2, ncol = 3)

  expect_equal(
    a %xn% (b ~ 1),
    n_prod(a, n = 1, mat = b)
  )
})

test_that("formula infix n-mode product preserves tensor class", {
  a <- tensor(array(1:3, dim = c(3, 1, 1)))
  b <- matrix(4:9, nrow = 2, ncol = 3)

  expect_s3_class(a %xn% (b ~ 1), "tensor")
})
