test_that("dtnorm computes large log determinants without underflow", {
  dims <- c(2, 50)
  x <- array(0, dim = dims)
  mu <- array(0, dim = dims)
  sigmas <- list(diag(2), diag(1e-10, 50))

  observed <- dtnorm(x, mu = mu, sigmas = sigmas, log = TRUE)
  full_covariance <- kronecker(sigmas[[2]], sigmas[[1]])
  expected <- mvtnorm::dmvnorm(
    c(x),
    mean = c(mu),
    sigma = full_covariance,
    log = TRUE
  )

  expect_true(is.finite(observed))
  expect_equal(observed, expected, tolerance = 1e-8)
})

test_that("tensor and matricized normal likelihoods agree", {
  set.seed(20260712)
  dims <- c(2, 3, 4)
  x <- array(rnorm(prod(dims)), dim = dims)
  mu <- array(rnorm(prod(dims)), dim = dims)
  sigmas <- lapply(dims, function(dimension) {
    z <- matrix(rnorm(dimension^2), dimension, dimension)
    crossprod(z) + diag(dimension)
  })

  tensor_loglik <- dtnorm(x, mu = mu, sigmas = sigmas, log = TRUE)

  for (mode in seq_along(dims)) {
    other_modes <- setdiff(seq_along(dims), mode)
    column_covariance <- Reduce(
      kronecker,
      rev(sigmas[other_modes])
    )

    matrix_loglik <- dtnorm(
      matricization(x, mode),
      mu = matricization(mu, mode),
      sigmas = list(sigmas[[mode]], column_covariance),
      log = TRUE
    )

    expect_equal(matrix_loglik, tensor_loglik, tolerance = 1e-8)
  }
})

test_that("dtnorm handles one tensor or all tensor draws", {
  dims <- c(2, 3)
  mu <- tensor(array(0, dim = dims))
  sigmas <- list(diag(2), diag(3))
  samples <- tensor(array(rnorm(5 * prod(dims)), dim = c(5, dims)), obs = 1)

  one <- dtnorm(pull_draw(samples, 1), mu = mu, sigmas = sigmas)
  all <- dtnorm(samples, mu = mu, sigmas = sigmas)
  expected <- vapply(
    seq_len(n_draws(samples)),
    function(i) dtnorm(pull_draw(samples, i), mu = mu, sigmas = sigmas),
    numeric(1)
  )

  expect_length(one, 1L)
  expect_length(all, n_draws(samples))
  expect_equal(all, expected)
})

test_that("dtnorm rejects list inputs", {
  x <- list(array(0, dim = c(2, 2)))

  expect_error(
    dtnorm(x, mu = array(0, dim = c(2, 2)), sigmas = list(diag(2), diag(2))),
    "lists are not supported"
  )
})
