test_that("univariate transform agrees with the Student-t CDF at zero skew", {
  values <- c(-2, 0, 1.5)
  nu <- 6
  draws <- tensor(values, dim = 1, n = length(values))

  transformed <- rosenblatt_tskewt(
    draws,
    mu = 0,
    skew = 0,
    sigmas = list(matrix(1)),
    nu = nu,
    rel.tol = 1e-9,
    abs.tol = 1e-11
  )

  expect_equal(
    as.numeric(unclass(transformed)),
    stats::pt(values, df = nu),
    tolerance = 1e-7
  )
})

test_that("modewise whitening gives the known symmetric-t Rosenblatt values", {
  nu <- 7
  dims <- c(2, 2)
  sigmas <- list(
    matrix(c(1.5, 0.3, 0.3, 0.9), 2, 2),
    matrix(c(0.8, -0.2, -0.2, 1.2), 2, 2)
  )
  lower_chols <- lapply(sigmas, function(sigma) t(chol(sigma)))

  standardized <- array(c(-0.8, 0.25, 1.1, -0.4), dim = dims)
  observed <- standardized
  for (mode in seq_along(lower_chols)) {
    observed <- n_prod(observed, mode, lower_chols[[mode]])
  }

  transformed <- rosenblatt_tskewt(
    tensor(observed),
    mu = array(0, dim = dims),
    skew = array(0, dim = dims),
    sigmas = sigmas,
    nu = nu,
    rel.tol = 1e-9,
    abs.tol = 1e-11
  )

  standardized_vector <- as.numeric(standardized)
  expected <- numeric(length(standardized_vector))
  cumulative_squared <- 0

  for (coordinate in seq_along(standardized_vector)) {
    previous <- coordinate - 1
    conditional_df <- nu + previous
    conditional_scale <- sqrt(
      (nu + cumulative_squared) / conditional_df
    )
    expected[coordinate] <- stats::pt(
      standardized_vector[coordinate] / conditional_scale,
      df = conditional_df
    )
    cumulative_squared <-
      cumulative_squared + standardized_vector[coordinate]^2
  }

  expect_equal(
    as.numeric(unclass(transformed)),
    expected,
    tolerance = 1e-7
  )
})

test_that("quadrature respects the skew-reflection identity", {
  lower <- tensormodels:::.tskewt_rosenblatt_conditional_cdf(
    value = 0.4,
    skew = 0.7,
    lambda = -4,
    chi = 5,
    psi = 0.6,
    rel.tol = 1e-9,
    abs.tol = 1e-11,
    subdivisions = 200L
  )
  reflected <- tensormodels:::.tskewt_rosenblatt_conditional_cdf(
    value = -0.4,
    skew = -0.7,
    lambda = -4,
    chi = 5,
    psi = 0.6,
    rel.tol = 1e-9,
    abs.tol = 1e-11,
    subdivisions = 200L
  )

  expect_equal(lower, 1 - reflected, tolerance = 1e-7)
})

test_that("transform preserves tensor dimensions and returns probabilities", {
  set.seed(42)
  mu <- array(0, dim = c(2, 2))
  skew <- array(c(0.3, -0.1, 0.2, 0.05), dim = c(2, 2))
  sigmas <- list(diag(2), diag(2))
  draws <- rtskewt(3, mu, sigmas, skew, nu = 8)

  transformed <- rosenblatt_tskewt(
    draws,
    mu = mu,
    skew = skew,
    sigmas = sigmas,
    nu = 8
  )

  expect_s3_class(transformed, "tensor")
  expect_equal(n_draws(transformed), n_draws(draws))
  expect_equal(draw_shape(transformed), draw_shape(draws))
  expect_true(all(is.finite(unclass(transformed))))
  expect_true(all(unclass(transformed) >= 0 & unclass(transformed) <= 1))
})

test_that("known-parameter skewed-t draws become approximately iid uniforms", {
  set.seed(20260730)
  mu <- array(c(-0.2, 0.3), dim = 2)
  skew <- array(c(0.6, -0.25), dim = 2)
  sigmas <- list(matrix(c(1.2, 0.35, 0.35, 0.9), 2, 2))
  draws <- rtskewt(
    300,
    mu = mu,
    skew = skew,
    sigmas = sigmas,
    nu = 7
  )

  transformed <- rosenblatt_tskewt(
    draws,
    mu = mu,
    skew = skew,
    sigmas = sigmas,
    nu = 7
  )
  uniforms <- matrix(
    unclass(transformed),
    nrow = n_draws(transformed)
  )

  expect_lt(max(abs(colMeans(uniforms) - 0.5)), 0.08)
  expect_lt(max(abs(apply(uniforms, 2, var) - 1 / 12)), 0.025)
  expect_lt(abs(stats::cor(uniforms)[1, 2]), 0.2)
})

test_that("transform validates its principal inputs", {
  draws <- tensor(array(0, dim = c(2, 2)))

  expect_error(
    rosenblatt_tskewt(
      array(0, dim = c(1, 1)),
      mu = array(0, dim = c(2, 2)),
      skew = array(0, dim = c(2, 2)),
      sigmas = list(diag(2), diag(2)),
      nu = 5
    ),
    "`draws` must be a `tensor`"
  )
  expect_error(
    rosenblatt_tskewt(
      draws,
      mu = array(0, dim = c(2, 2)),
      skew = array(0, dim = c(2, 2)),
      sigmas = list(diag(2), diag(2)),
      nu = 0
    ),
    "`parameters\\$nu` must be one positive"
  )
})
