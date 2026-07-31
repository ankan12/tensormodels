test_that("score-adjusted GOF test returns a finite htest for tensor normals", {
  set.seed(20260731)
  mu <- array(c(-0.25, 0.5), dim = 2)
  sigmas <- list(matrix(c(1.2, 0.2, 0.2, 0.9), 2, 2))
  draws <- rtnorm(120, mu = mu, sigmas = sigmas)

  result <- tensor_score_gof(
    draws,
    model = "normal",
    marginal_degree = 2,
    include_interactions = TRUE,
    max_basis = 5
  )

  expect_s3_class(result, "htest")
  expect_s3_class(result, "tensor_score_gof_test")
  expect_true(is.finite(unname(result$statistic)))
  expect_true(result$p.value >= 0 && result$p.value <= 1)
  expect_equal(result$score_dimension, 5L)
  expect_equal(result$information_rank, 5L)
  expect_length(result$basis_names, 5L)
  expect_equal(dim(unclass(result$uniforms)), c(120, 2))
})

test_that("smooth basis contains marginal and dependence functions", {
  uniforms <- matrix(
    c(
      0.1, 0.3, 0.7, 0.9,
      0.2, 0.4, 0.6, 0.8
    ),
    nrow = 4,
    ncol = 2
  )

  basis <- make_basis(
    uniforms,
    marginal_degree = 2,
    include_interactions = TRUE,
    max_basis = 5
  )

  expect_equal(
    colnames(basis),
    c(
      "L1(U1)", "L1(U2)", "L2(U1)", "L2(U2)",
      "L1(U1)*L1(U2)"
    )
  )
  expect_equal(dim(basis), c(4, 5))
  expect_equal(
    basis[, 5],
    basis[, 1] * basis[, 2]
  )
})

test_that("make_basis accepts tensors, arrays, and vectors", {
  values <- array(
    seq(0.05, 0.95, length.out = 24),
    dim = c(4, 2, 3),
    dimnames = list(paste0("draw", 1:4), NULL, NULL)
  )

  array_basis <- make_basis(
    values,
    marginal_degree = 1,
    include_interactions = FALSE
  )
  tensor_basis <- make_basis(
    tensor(values, obs = 1),
    marginal_degree = 1,
    include_interactions = FALSE
  )
  vector_basis <- make_basis(
    c(0.2, 0.4, 0.6),
    marginal_degree = 2,
    include_interactions = TRUE
  )

  expect_equal(dim(array_basis), c(4, 6))
  expect_equal(tensor_basis, array_basis)
  expect_equal(dim(vector_basis), c(3, 2))
  expect_equal(
    colnames(vector_basis),
    c("L1(U1)", "L2(U1)")
  )
})

test_that("make_basis validates uniform values", {
  expect_error(
    make_basis(c(0.2, 1.1)),
    "between zero and one"
  )
  expect_error(
    make_basis("not numeric"),
    "numeric vector"
  )
})

test_that("finite-difference scores support every fitted distribution", {
  set.seed(20260801)
  mu <- array(0, dim = 1)
  skew <- array(0.3, dim = 1)
  sigmas <- list(matrix(1.2))

  specifications <- list(
    normal = list(
      draws = rtnorm(8, mu = mu, sigmas = sigmas),
      fit = list(mu = mu, sigmas = sigmas),
      score_dimension = 2L
    ),
    skewt = list(
      draws = rtskewt(
        8, mu = mu, sigmas = sigmas, skew = skew, nu = 7
      ),
      fit = list(mu = mu, skew = skew, sigmas = sigmas, nu = 7),
      score_dimension = 4L
    ),
    vargamma = list(
      draws = rtvargamma(
        8, mu = mu, sigmas = sigmas, skew = skew, scale = 4
      ),
      fit = list(mu = mu, skew = skew, sigmas = sigmas, gamma = 4),
      score_dimension = 4L
    ),
    invgauss = list(
      draws = rtinvgauss(
        8, mu = mu, sigmas = sigmas, skew = skew, kappa = 2
      ),
      fit = list(mu = mu, skew = skew, sigmas = sigmas, kappa = 2),
      score_dimension = 4L
    ),
    genhyper = list(
      draws = rtgenhyper(
        8, mu = mu, sigmas = sigmas, skew = skew,
        lambda = 1, omega = 4
      ),
      fit = list(
        mu = mu, skew = skew, sigmas = sigmas,
        lambda = 1, omega = 4
      ),
      score_dimension = 5L
    )
  )

  for (model in names(specifications)) {
    specification <- specifications[[model]]
    parameterization <-
      tensormodels:::.tensor_score_gof_parameterization(
        specification$fit,
        model,
        dims = 1L
      )
    scores <- tensormodels:::.tensor_score_gof_scores(
      specification$draws,
      model,
      parameterization,
      fd_step = 1e-5
    )

    expect_equal(
      ncol(scores),
      specification$score_dimension,
      info = model
    )
    expect_true(all(is.finite(scores)), info = model)
  }
})

test_that("score-adjusted GOF validates unsupported configurations", {
  draws <- rtnorm(5, mu = array(0, dim = 1))

  expect_error(
    tensor_score_gof(draws, model = "normal", marginal_degree = 0),
    "positive integer"
  )
  expect_error(
    tensor_score_gof(draws, model = "normal", restrict = 1),
    "not currently supported"
  )
  expect_error(
    tensor_score_gof(unclass(draws), model = "normal"),
    "must be a `tensor`"
  )
})
