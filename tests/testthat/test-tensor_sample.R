make_tensor_samples <- function(n, shape) {
  input <- array(seq_len(n * prod(shape)), dim = c(shape, n))
  tensor_samples(input, obs = length(shape) + 1L)
}

test_that("tensor_samples constructs an array-backed sample", {
  input <- array(1:24, dim = c(3, 4, 2))
  x <- tensor_samples(input, obs = 3)

  expect_s3_class(x, "tensor_samples")
  expect_true(is.array(x))
  expect_equal(dim(x), c(2L, 3L, 4L))
  expect_equal(vctrs::vec_size(x), 2L)
  expect_equal(n_draws(x), 2L)
  expect_equal(draw_shape(x), c(3L, 4L))
})

test_that("tensor_samples requires an array and validates the observation axis", {
  expect_error(tensor_samples(1:12), "at least two dimensions")
  expect_error(tensor_samples(array(1:4, dim = 4)), "at least two dimensions")
  expect_error(tensor_samples(array(1:12, c(2, 2, 3)), obs = 4), "past the end")
  expect_error(tensor_samples(array(1:12, c(2, 2, 3)), obs = 1.5), "loss of precision")
})

test_that("tensor_samples moves the observation axis to the front", {
  first <- array(1:12, dim = c(3, 4))
  second <- array(13:24, dim = c(3, 4))
  input <- simplify2array(list(first, second))

  x <- tensor_samples(input, obs = 3)

  expect_equal(dim(x), c(2L, 3L, 4L))
  expect_equal(unclass(pull_draw(x, 1)), first)
  expect_equal(unclass(pull_draw(x, 2)), second)
})

test_that("tensor_samples accepts a named observation axis", {
  input <- array(1:24, dim = c(row = 3, column = 4, person = 2))
  x <- tensor_samples(input, obs = "person")

  expect_equal(names(dim(x)), c("person", "row", "column"))
  expect_equal(dim(x), c(person = 2L, row = 3L, column = 4L))
})

test_that("tensor_samples has an informative compact print method", {
  x <- make_tensor_samples(2, c(3, 4))

  expect_output(print(x), "<tensor_samples\\[2\\]>")
  expect_output(print(x), "2 observations of shape 3 × 4", fixed = TRUE)
  expect_output(print(x), "Preview: 2 observations × 6 indexed entries", fixed = TRUE)
  expect_output(print(x), "[1,1] [2,1] [3,1] [1,2] [2,2] [3,2]", fixed = TRUE)
  expect_output(print(x), "6 more entries per observation", fixed = TRUE)
  expect_identical(print(x), x)
})

test_that("tensor_samples print preview is bounded and configurable", {
  x <- make_tensor_samples(10, c(2, 3))

  expect_output(
    print(x, n = 2, entries = 3),
    "Preview: 2 observations × 3 indexed entries",
    fixed = TRUE
  )
  expect_output(
    print(x, n = 2, entries = 3),
    "8 more observations and 3 more entries per observation",
    fixed = TRUE
  )
  expect_output(print(x, n = 2, entries = 3), "obs 1")
  expect_output(print(x, n = 2, entries = 3), "obs 2")
})

test_that("tensor_samples print preview uses observation names", {
  x <- make_tensor_samples(2, c(2, 2))
  dimnames(x) <- list(c("person-a", "person-b"), NULL, NULL)

  expect_output(print(x), "person-a")
  expect_output(print(x), "person-b")
})

test_that("pull_draw returns one tensor without the observation axis", {
  x <- make_tensor_samples(2, c(3, 4))
  first <- pull_draw(x, 1)

  expect_false(inherits(first, "tensor_samples"))
  expect_s3_class(first, "tensor")
  expect_true(is.array(first))
  expect_equal(dim(first), c(3L, 4L))
  expect_equal(unclass(first), unclass(x)[1, , ])
})

test_that("pull_draw preserves tensor mode dimnames", {
  x <- make_tensor_samples(2, c(2, 3))
  dimnames(x) <- list(
    draw = c("first", "second"),
    row = c("r1", "r2"),
    column = c("c1", "c2", "c3")
  )

  second <- pull_draw(x, "second")

  expect_equal(dimnames(second), dimnames(x)[-1L])
  expect_equal(unclass(second), unclass(x)[2, , ])
})

test_that("slice_draws retains the tensor samples class and shape", {
  x <- make_tensor_samples(4, c(2, 3))
  selected <- slice_draws(x, c(4, 2))

  expect_s3_class(selected, "tensor_samples")
  expect_equal(dim(selected), c(2L, 2L, 3L))
  expect_equal(n_draws(selected), 2L)
  expect_equal(draw_shape(selected), c(2L, 3L))
})

test_that("vctrs combines compatible tensor samples along the draw axis", {
  x <- make_tensor_samples(2, c(2, 3))
  y <- make_tensor_samples(3, c(2, 3))

  combined <- vctrs::vec_c(x, y)

  expect_s3_class(combined, "tensor_samples")
  expect_equal(dim(combined), c(5L, 2L, 3L))
  expect_equal(n_draws(combined), 5L)
})
