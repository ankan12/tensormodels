test_that("tensor stores one draw with an explicit observation axis", {
  input <- array(1:24, dim = c(2, 3, 4))
  x <- tensor(input)

  expect_s3_class(x, "tensor")
  expect_true(is.array(x))
  expect_equal(dim(x), c(2L, 3L, 4L))
  expect_equal(n_draws(x), 1L)
  expect_equal(draw_shape(x), c(2L, 3L, 4L))
  expect_equal(.tensor_single_draw_array(pull_draw(x, 1)), input)
})

test_that("tensor accepts array-style data and dimensions", {
  x <- tensor(data = 1:24, dim = c(4, 3, 2))

  expect_equal(n_draws(x), 1L)
  expect_equal(draw_shape(x), c(4L, 3L, 2L))
  expect_equal(.tensor_single_draw_array(pull_draw(x, 1)),
               array(1:24, dim = c(4, 3, 2)))
  expect_warning(
    tensor(data = 1:12, dim = c(4, 3, 2)),
    "values will be recycled"
  )
  expect_equal(x[3, 1, 1], 3)
  x[3, 1, 1] <- 99
  expect_equal(x[3, 1, 1], 99)
})

test_that("tensor moves an explicit observation axis", {
  input <- array(1:48, dim = c(2, 3, 4, 2))
  x <- tensor(input, obs = 4)

  expect_equal(dim(x), c(2L, 3L, 4L))
  expect_equal(n_draws(x), 2L)
  expect_equal(draw_shape(x), c(2L, 3L, 4L))
})

test_that("tensor arithmetic broadcasts a singleton draw", {
  draws <- tensor(array(1:48, dim = c(2, 3, 4, 2)), obs = 4)
  offset <- tensor(array(rep(1, 24), dim = c(2, 3, 4)))
  result <- draws + offset

  expect_s3_class(result, "tensor")
  expect_equal(dim(result), dim(draws))
  expect_equal(.tensor_single_draw_array(pull_draw(result, 1)),
               .tensor_single_draw_array(pull_draw(draws, 1)) + 1)
  expect_equal(.tensor_single_draw_array(pull_draw(result, 2)),
               .tensor_single_draw_array(pull_draw(draws, 2)) + 1)
})

test_that("tensor arithmetic rejects incompatible shapes and draw counts", {
  x <- tensor(array(1:4, dim = c(2, 2)))
  y <- tensor(array(1:6, dim = c(2, 3)))
  z <- tensor(array(1:12, dim = c(2, 2, 3)), obs = 3)

  expect_error(x + y, "draw shapes must agree")
  expect_error(z + tensor(array(1:8, dim = c(2, 2, 2)), obs = 3),
               "draw counts must agree")
})

test_that("as_tensor uses the tensor container", {
  x <- as_tensor(array(1:6, dim = c(2, 3)))
  expect_s3_class(x, "tensor")
  expect_equal(n_draws(x), 1L)
  expect_equal(draw_shape(x), c(2L, 3L))
})
