test_that("analytic tensor scores agree with likelihood derivatives", {
  set.seed(20260802)
  dims <- c(2L, 2L)
  mu <- array(c(-0.2, 0.1, 0.4, -0.1), dim = dims)
  skew <- array(c(0.2, -0.1, 0.1, 0.3), dim = dims)
  sigmas <- list(
    matrix(c(1, 0.15, 0.15, 1), 2, 2),
    matrix(c(1.3, 0.1, 0.1, 0.8), 2, 2)
  )

  specifications <- list(
    normal = list(
      draws = rtnorm(5, mu = mu, sigmas = sigmas),
      fit = list(mu = mu, sigmas = sigmas),
      tolerance = 1e-5
    ),
    skewt = list(
      draws = rtskewt(
        5, mu = mu, sigmas = sigmas, skew = skew, nu = 8
      ),
      fit = list(mu = mu, skew = skew, sigmas = sigmas, nu = 8),
      tolerance = 2e-4
    ),
    vargamma = list(
      draws = rtvargamma(
        5, mu = mu, sigmas = sigmas, skew = skew, scale = 5
      ),
      fit = list(mu = mu, skew = skew, sigmas = sigmas, gamma = 5),
      tolerance = 2e-3
    ),
    invgauss = list(
      draws = rtinvgauss(
        5, mu = mu, sigmas = sigmas, skew = skew, kappa = 2
      ),
      fit = list(mu = mu, skew = skew, sigmas = sigmas, kappa = 2),
      tolerance = 5e-3
    ),
    genhyper = list(
      draws = rtgenhyper(
        5, mu = mu, sigmas = sigmas, skew = skew,
        lambda = 1, omega = 4
      ),
      fit = list(
        mu = mu, skew = skew, sigmas = sigmas,
        lambda = 1, omega = 4
      ),
      tolerance = 5e-3
    )
  )

  for (model in names(specifications)) {
    specification <- specifications[[model]]
    parameterization <-
      tensortools:::.tensor_score_gof_parameterization(
        specification$fit,
        model,
        dims
      )
    analytic <- tensor_scores(
      specification$draws,
      model = model,
      fit = specification$fit
    )
    numerical <- tensortools:::.tensor_score_gof_scores(
      specification$draws,
      model,
      parameterization,
      fd_step = 1e-5
    )

    expect_equal(
      as.numeric(analytic),
      as.numeric(numerical),
      tolerance = specification$tolerance,
      info = model
    )
    expect_identical(colnames(analytic), parameterization$names)
    expect_identical(attr(analytic, "model"), model)
  }
})

test_that("tensor_scores handles zero skew in the skew-t model", {
  set.seed(20260803)
  mu <- array(0, dim = c(2L, 2L))
  skew <- array(0, dim = c(2L, 2L))
  sigmas <- list(diag(2), diag(2))
  draws <- rtskewt(
    10,
    mu = mu,
    sigmas = sigmas,
    skew = skew,
    nu = 7
  )

  scores <- tensor_scores(
    draws,
    model = "skewt",
    fit = list(mu = mu, skew = skew, sigmas = sigmas, nu = 7)
  )

  expect_equal(dim(scores), c(10L, 14L))
  expect_true(all(is.finite(scores)))
})

test_that("fast score-adjusted GOF returns statistic and p-value", {
  set.seed(20260804)
  mu <- array(c(-0.25, 0.5), dim = 2L)
  sigmas <- list(matrix(c(1.2, 0.2, 0.2, 0.9), 2, 2))
  draws <- rtnorm(120, mu = mu, sigmas = sigmas)

  result <- tensor_gof_scores(
    draws,
    model = "normal",
    marginal_degree = 2,
    include_interactions = TRUE,
    max_basis = 5,
    B = 99,
    seed = 1
  )

  expect_s3_class(result, "htest")
  expect_s3_class(result, "tensor_gof_scores_test")
  expect_true(is.finite(unname(result$statistic)))
  expect_true(result$p.value >= 0 && result$p.value <= 1)
  expect_equal(result$score_dimension, 5L)
  expect_equal(dim(result$scores), c(120L, 5L))
  expect_identical(result$calibration, "multiplier")
  expect_length(result$multiplier_statistics, 99L)
  expect_equal(
    unname(result$statistic),
    n_draws(draws) * sum(result$basis_means^2)
  )
  expect_true(is.finite(unname(result$asymptotic.statistic)))
  expect_true(
    result$asymptotic.p.value >= 0 &&
      result$asymptotic.p.value <= 1
  )
})

test_that("weighted multiplier calibration is reproducible", {
  set.seed(20260805)
  draws <- rtnorm(
    60,
    mu = array(0, dim = 2L),
    sigmas = list(diag(2))
  )
  fit <- tensor_mle(draws, model = "normal")

  first <- tensor_gof_scores(
    draws,
    model = "normal",
    fit = fit,
    max_basis = 4,
    B = 49,
    seed = 17,
    multiplier_chunk_size = 11
  )
  second <- tensor_gof_scores(
    draws,
    model = "normal",
    fit = fit,
    max_basis = 4,
    B = 49,
    seed = 17,
    multiplier_chunk_size = 17
  )

  expect_equal(first$p.value, second$p.value)
  expect_equal(
    first$multiplier_statistics,
    second$multiplier_statistics
  )
})

test_that("new score functions validate their inputs", {
  draws <- rtnorm(5, mu = array(0, dim = 1L))

  expect_error(
    tensor_scores(unclass(draws), model = "normal"),
    "must be a `tensor`"
  )
  expect_error(
    tensor_scores(draws, model = "normal", restrict = 1),
    "not currently supported"
  )
  expect_error(
    tensor_gof_scores(draws, model = "normal", restrict = 1),
    "not currently supported"
  )
  expect_error(
    tensor_gof_scores(draws, model = "normal", B = 0),
    "positive integer"
  )
})
