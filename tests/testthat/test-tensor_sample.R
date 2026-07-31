make_tensor <- function(n, shape) {
  input <- array(seq_len(n * prod(shape)), dim = c(shape, n))
  tensor(input, obs = length(shape) + 1L)
}

test_that("tensor prints a compact draw preview", {
  x <- make_tensor(2, c(3, 4))
  expect_output(print(x), "<tensor\\[2\\]>")
  expect_output(print(x), "2 draws of shape 3 × 4", fixed = TRUE)
  expect_output(print(x), "Preview: 2 draws × 6 indexed entries", fixed = TRUE)
  expect_output(print(x), "pull_draw", fixed = TRUE, negate = TRUE)
})

test_that("tensor subsetting and pull_draw preserve tensor structure", {
  x <- make_tensor(4, c(2, 3))
  selected <- slice_draws(x, c(4, 2))
  first <- pull_draw(x, 1)

  expect_s3_class(selected, "tensor")
  expect_equal(dim(selected), c(2L, 3L))
  expect_equal(n_draws(selected), 2L)
  expect_equal(dim(first), c(2L, 3L))
  expect_equal(n_draws(first), 1L)
  expect_equal(unclass(pull_draw(first, 1)), unclass(pull_draw(x, 1)))
})
