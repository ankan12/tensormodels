test_that("2x1 prod with matrix × matrix ", {
  A <- matrix(c(1, 2, 3, 4), nrow = 2)

  b <- matrix(c(5, 6), nrow = 2)

  out <- nm_prod(A, b, 2, 1)

  expect_equal(dim(out), c(2,1))
  expect_equal(out, matrix(c(23, 34), nrow = 2))
})

test_that("1x1 prod with matrix × matrix ", {
  A <- matrix(c(1, 2, 3, 4), nrow = 2)

  b <- matrix(c(5, 6), nrow = 2)

  out <- nm_prod(A, b, 1, 1)

  expect_equal(dim(out), c(2,1))
  expect_equal(out, matrix(c(17, 39), nrow = 2))
})
