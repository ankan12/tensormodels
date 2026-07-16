test_that("tensor wraps an array without changing its values or shape", {
  input <- array(1:24, dim = c(2, 3, 4))
  x <- tensor(input)

  expect_s3_class(x, "tensor")
  expect_true(is.array(x))
  expect_equal(dim(x), dim(input))
  expect_equal(unclass(x), input)
  expect_error(tensor(1:24), "Supply `dim`")
})

test_that("tensor constructs an array from data and dimensions", {
  input <- array(1:24, dim = c(2, 3, 4))
  x <- tensor(1:24, dim = c(2, 3, 4))

  expect_s3_class(x, "tensor")
  expect_equal(unclass(x), input)
  expect_equal(tensor(x = 1:24, dim = c(2, 3, 4)), x)
})

test_that("tensor validates constructor dimensions", {
  expect_error(
    tensor(1:5, dim = c(2, 3)),
    "length of `x` must equal"
  )
  expect_error(
    tensor(1:6, dim = c(2, 0, 3)),
    "positive integers"
  )
})

test_that("tensor arithmetic is elementwise and preserves the class", {
  x <- tensor(array(1:4, dim = c(2, 2)))
  y <- tensor(array(5:8, dim = c(2, 2)))

  expect_s3_class(x + y, "tensor")
  expect_equal(unclass(x + y), unclass(x) + unclass(y))
  expect_equal(unclass(x - y), unclass(x) - unclass(y))
  expect_equal(unclass(x * y), unclass(x) * unclass(y))
  expect_equal(unclass(x / 2), unclass(x) / 2)
  expect_equal(unclass(2 * x), 2 * unclass(x))
  expect_equal(unclass(-x), -unclass(x))
})

test_that("tensor arithmetic rejects incompatible tensor shapes", {
  x <- tensor(array(1:4, dim = c(2, 2)))
  y <- tensor(array(1:6, dim = c(2, 3)))

  expect_error(x + y, "Tensor shapes must agree: 2 x 2 and 2 x 3", fixed = TRUE)
})

test_that("tensor comparisons return ordinary logical arrays", {
  x <- tensor(array(1:4, dim = c(2, 2)))
  out <- x > 2

  expect_false(inherits(out, "tensor"))
  expect_true(is.array(out))
  expect_type(out, "logical")
  expect_equal(out, unclass(x) > 2)
})

test_that("as_tensor uses the general tensor constructor", {
  input <- array(1:6, dim = c(2, 3))
  x <- as_tensor(input)

  expect_s3_class(x, "tensor")
  expect_equal(unclass(x), input)
})
