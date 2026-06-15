test_that("plot_density returns a ggplot for normal density slices", {
  mu <- matrix(0, nrow = 2, ncol = 2)
  sigmas <- list(diag(2), diag(2))

  p <- plot_density(
    x12 = c(-1, 0, 1),
    x22 = c(-1, 0, 1),
    x11 = seq(-2, 2, length.out = 8),
    x21 = seq(-2, 2, length.out = 8),
    mu = mu,
    sigmas = sigmas,
    model = "normal"
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_density returns a ggplot for skewed density slices", {
  mu <- matrix(0, nrow = 2, ncol = 2)
  sigmas <- list(diag(2), diag(2))
  skew <- matrix(0.1, nrow = 2, ncol = 2)

  p <- plot_density(
    x12 = c(-1, 0, 1),
    x22 = c(-1, 0, 1),
    x11 = seq(-2, 2, length.out = 8),
    x21 = seq(-2, 2, length.out = 8),
    mu = mu,
    sigmas = sigmas,
    model = "skewt",
    skew = skew,
    nu = 20
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_density checks required skewed model parameters", {
  mu <- matrix(0, nrow = 2, ncol = 2)
  sigmas <- list(diag(2), diag(2))

  expect_error(
    plot_density(
      x12 = c(-1, 0, 1),
      x22 = c(-1, 0, 1),
      x11 = seq(-2, 2, length.out = 8),
      x21 = seq(-2, 2, length.out = 8),
      mu = mu,
      sigmas = sigmas,
      model = "skewt",
      nu = 20
    ),
    "skew must be supplied"
  )
})
