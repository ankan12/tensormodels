test_that("normal transform recovers exact modewise Gaussian innovations", {
  dims <- c(2, 2)
  mu <- array(c(-0.5, 0, 0.5, 1), dim = dims)
  sigmas <- list(
    matrix(c(1.4, 0.25, 0.25, 0.8), 2, 2),
    matrix(c(0.9, -0.15, -0.15, 1.2), 2, 2)
  )
  lower_chols <- lapply(sigmas, function(sigma) t(chol(sigma)))
  innovations <- array(c(-1.2, 0.3, 0.8, -0.4), dim = dims)

  observed <- innovations
  for (mode in seq_along(lower_chols)) {
    observed <- n_prod(observed, mode, lower_chols[[mode]])
  }
  observed <- observed + mu

  transformed <- tensor_rosenblatt(
    tensor(observed),
    model = "normal",
    parameters = list(mu = mu, sigmas = sigmas)
  )

  expect_equal(
    as.numeric(unclass(transformed)),
    stats::pnorm(as.numeric(innovations)),
    tolerance = 1e-10
  )
})

test_that("all supported known-parameter models give approximately iid uniforms", {
  mu <- array(c(-0.2, 0.3), dim = 2)
  skew <- array(c(0.45, -0.2), dim = 2)
  sigmas <- list(matrix(c(1.2, 0.3, 0.3, 0.9), 2, 2))

  specifications <- list(
    skewt = list(
      simulate = function() {
        rtskewt(300, mu, sigmas, skew, nu = 7)
      },
      parameters = list(
        mu = mu, skew = skew, sigmas = sigmas, nu = 7
      )
    ),
    vargamma = list(
      simulate = function() {
        rtvargamma(300, mu, sigmas, skew, scale = 4)
      },
      parameters = list(
        mu = mu, skew = skew, sigmas = sigmas, gamma = 4
      )
    ),
    invgauss = list(
      simulate = function() {
        rtinvgauss(300, mu, sigmas, skew, kappa = 2)
      },
      parameters = list(
        mu = mu, skew = skew, sigmas = sigmas, kappa = 2
      )
    ),
    genhyper = list(
      simulate = function() {
        rtgenhyper(
          300, mu = mu, skew = skew, sigmas = sigmas,
          lambda = 1, omega = 4
        )
      },
      parameters = list(
        mu = mu, skew = skew, sigmas = sigmas,
        lambda = 1, omega = 4
      )
    )
  )

  for (model_index in seq_along(specifications)) {
    model <- names(specifications)[model_index]
    specification <- specifications[[model]]
    set.seed(20260730 + model_index)
    draws <- specification$simulate()

    transformed <- tensor_rosenblatt(
      draws,
      model = model,
      parameters = specification$parameters
    )
    uniforms <- matrix(
      unclass(transformed),
      nrow = n_draws(transformed)
    )

    expect_lt(
      max(abs(colMeans(uniforms) - 0.5)),
      0.08
    )
    expect_lt(
      max(abs(apply(uniforms, 2, var) - 1 / 12)),
      0.025
    )
    expect_lt(
      abs(stats::cor(uniforms)[1, 2]),
      0.2
    )
  }
})

test_that("variance-gamma accepts sampler and MLE parameter names", {
  set.seed(17)
  mu <- array(c(0, 0), dim = 2)
  skew <- array(c(0.2, -0.1), dim = 2)
  sigmas <- list(diag(2))
  draws <- rtvargamma(5, mu, sigmas, skew, scale = 3)

  with_gamma <- tensor_rosenblatt(
    draws,
    model = "vargamma",
    parameters = list(
      mu = mu, skew = skew, sigmas = sigmas, gamma = 3
    )
  )
  with_scale <- tensor_rosenblatt(
    draws,
    model = "vargamma",
    parameters = list(
      mu = mu, skew = skew, sigmas = sigmas, scale = 3
    )
  )

  expect_equal(unclass(with_gamma), unclass(with_scale))
})

test_that("specialized skew-t wrapper delegates to unified transform", {
  set.seed(21)
  mu <- array(c(0, 0), dim = 2)
  skew <- array(c(0.3, -0.2), dim = 2)
  sigmas <- list(diag(2))
  draws <- rtskewt(5, mu, sigmas, skew, nu = 6)

  unified <- tensor_rosenblatt(
    draws,
    model = "skewt",
    parameters = list(
      mu = mu, skew = skew, sigmas = sigmas, nu = 6
    )
  )
  specialized <- rosenblatt_tskewt(
    draws,
    mu = mu,
    skew = skew,
    sigmas = sigmas,
    nu = 6
  )

  expect_equal(unclass(unified), unclass(specialized))
})

test_that("unified transform validates model-specific parameter lists", {
  draws <- tensor(array(0, dim = c(2, 2)))
  common <- list(
    mu = array(0, dim = c(2, 2)),
    skew = array(0, dim = c(2, 2)),
    sigmas = list(diag(2), diag(2))
  )

  expect_error(
    tensor_rosenblatt(
      draws,
      model = "skewt",
      parameters = common
    ),
    "must contain `nu`"
  )
  expect_error(
    tensor_rosenblatt(
      draws,
      model = "vargamma",
      parameters = c(common, list(gamma = 2, scale = 3))
    ),
    "must agree"
  )
  expect_error(
    tensor_rosenblatt(
      draws,
      model = "genhyper",
      parameters = c(common, list(lambda = NA_real_, omega = 2))
    ),
    "`parameters\\$lambda`"
  )
})
