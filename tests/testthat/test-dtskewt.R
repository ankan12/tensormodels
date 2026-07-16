test_that("density functions check parameter dimensions", {
  x <- array(1:6, dim = c(2, 3))
  mu <- array(0, dim = c(2, 3))
  skew <- array(1, dim = c(2, 3))
  sigmas <- list(diag(2), diag(3))
  skew_densities <- list(
    dtinvgauss = function(...) dtinvgauss(..., kappa = 2),
    dtskewt = function(...) dtskewt(..., nu = 4),
    dtvargamma = function(...) dtvargamma(..., scale = 2),
    dtgenhyper = function(...) dtgenhyper(..., lambda = 2, omega = 2)
  )

  expect_error(
    dtnorm(x, mu = array(0, dim = c(3, 2)), sigmas = sigmas),
    "mu must have the same dimensions as x",
    fixed = TRUE
  )

  expect_error(
    dtnorm(x, mu = mu, sigmas = list(diag(2))),
    "sigmas must be a list with one covariance matrix for each dimension",
    fixed = TRUE
  )

  expect_error(
    dtnorm(x, mu = mu, sigmas = list(diag(2), diag(2))),
    "sigmas[[2]] must be a 3 x 3 covariance matrix",
    fixed = TRUE
  )

  for(density in skew_densities) {
    expect_error(
      density(x, mu = array(0, dim = c(3, 2)), skew = skew,
              sigmas = sigmas),
      "mu must have the same dimensions as x",
      fixed = TRUE
    )

    expect_error(
      density(x, mu = mu, skew = array(1, dim = c(3, 2)),
              sigmas = sigmas),
      "skew must have the same dimensions as x",
      fixed = TRUE
    )

    expect_error(
      density(x, mu = mu, skew = skew, sigmas = list(diag(2))),
      "sigmas must be a list with one covariance matrix for each dimension",
      fixed = TRUE
    )

    expect_error(
      density(x, mu = mu, skew = skew, sigmas = list(diag(2), diag(2))),
      "sigmas[[2]] must be a 3 x 3 covariance matrix",
      fixed = TRUE
    )
  }
})

test_that("density functions accept vector-valued inputs", {
  x <- c(1, 2)
  mu <- c(0, 0)
  skew <- c(1, 1)
  sigmas <- list(diag(2))

  expect_type(dtnorm(x, mu = mu, sigmas = sigmas), "double")
  expect_type(dtinvgauss(x, mu = mu, skew = skew, sigmas = sigmas,
                         kappa = 2), "double")
  expect_type(dtskewt(x, mu = mu, skew = skew, sigmas = sigmas,
                      nu = 4), "double")
  expect_type(dtvargamma(x, mu = mu, skew = skew, sigmas = sigmas,
                         scale = 2), "double")
  expect_type(dtgenhyper(x, mu = mu, skew = skew, sigmas = sigmas,
                         lambda = 2, omega = 2), "double")
})

test_that("skewed density functions warn and regularize singular covariance log determinants", {
  x <- c(1, 2)
  mu <- c(0, 0)
  skew <- c(1, 1)
  sigmas <- list(matrix(c(1, 0, 0, 0), nrow = 2))
  skew_densities <- list(
    dtinvgauss = function(...) dtinvgauss(..., kappa = 2, log = TRUE),
    dtskewt = function(...) dtskewt(..., nu = 4, log = TRUE),
    dtvargamma = function(...) dtvargamma(..., scale = 2, log = TRUE),
    dtgenhyper = function(...) dtgenhyper(..., lambda = 2, omega = 2, log = TRUE)
  )

  for (density in skew_densities) {
    expect_warning(
      value <- density(x, mu = mu, skew = skew, sigmas = sigmas),
      "ridge regularization",
      fixed = TRUE
    )
    expect_true(is.finite(value))
  }
})

test_that("random generation functions check parameter dimensions", {
  mu <- array(0, dim = c(2, 3))
  skew <- array(1, dim = c(2, 3))
  sigmas <- list(diag(2), diag(3))
  skew_generators <- list(
    rtinvgauss = function(...) rtinvgauss(..., kappa = 2),
    rtskewt = function(...) rtskewt(..., nu = 4),
    rtvargamma = function(...) rtvargamma(..., scale = 2),
    rtgenhyper = function(...) rtgenhyper(..., lambda = 2, omega = 2)
  )

  expect_error(
    rtnorm(1, mu = mu, sigmas = list(diag(2))),
    "sigmas must be a list with one covariance matrix for each dimension",
    fixed = TRUE
  )

  expect_error(
    rtnorm(1, mu = mu, sigmas = list(diag(2), diag(2))),
    "sigmas[[2]] must be a 3 x 3 covariance matrix",
    fixed = TRUE
  )

  for(generator in skew_generators) {
    expect_error(
      generator(1, mu = mu, skew = array(1, dim = c(3, 2)),
                sigmas = sigmas),
      "skew must have the same dimensions as mu",
      fixed = TRUE
    )

    expect_error(
      generator(1, mu = mu, skew = skew, sigmas = list(diag(2))),
      "sigmas must be a list with one covariance matrix for each dimension",
      fixed = TRUE
    )

    expect_error(
      generator(1, mu = mu, skew = skew, sigmas = list(diag(2), diag(2))),
      "sigmas[[2]] must be a 3 x 3 covariance matrix",
      fixed = TRUE
    )
  }
})

test_that("random generation functions accept matching tensor parameters", {
  mu <- array(0, dim = c(2, 3))
  skew <- array(1, dim = c(2, 3))
  sigmas <- list(diag(2), diag(3))

  expect_equal(dim(rtnorm(1, mu = mu, sigmas = sigmas)[[1]]), c(2, 3))
  expect_equal(dim(rtinvgauss(1, mu = mu, skew = skew, sigmas = sigmas,
                              kappa = 2)[[1]]), c(2, 3))
  expect_equal(dim(rtskewt(1, mu = mu, skew = skew, sigmas = sigmas,
                           nu = 4)[[1]]), c(2, 3))
  expect_equal(dim(rtvargamma(1, mu = mu, skew = skew, sigmas = sigmas,
                              scale = 2)[[1]]), c(2, 3))
  expect_equal(dim(rtgenhyper(1, mu = mu, skew = skew, sigmas = sigmas,
                              lambda = 2, omega = 2)[[1]]), c(2, 3))
})
