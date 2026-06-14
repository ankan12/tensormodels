.eval_density_input <- function(x, eval_one, log) {
  if (is.list(x)) {
    return(vapply(x, eval_one, numeric(1), log = log))
  }

  eval_one(x, log = log)
}

.density_all_det <- function(sigmas, dims, zero_floor = FALSE) {
  n_star <- prod(dims)
  all_det <- 0

  for (k in seq_along(sigmas)) {
    sigd <- sigmas[[k]]
    det_sigd <- det(sigd)

    if (zero_floor && det_sigd == 0) det_sigd <- 1e-8

    all_det <- all_det +
      (n_star / (2 * nrow(sigd))) * log(det_sigd)
  }

  all_det
}

.make_dtnorm_evaluator <- function(mu, sigmas, log = FALSE) {
  default_log <- log
  dims <- .tensor_dims(mu)
  n_star <- prod(dims)

  .validate_same_dims(mu, dims, "mu")
  sigmas <- .prepare_sigmas(sigmas, dims)

  mu <- array(mu, dim = dims)
  inv_sigmas <- lapply(sigmas, invert_safe)
  all_det <- .density_all_det(sigmas, dims, zero_floor = TRUE)
  const <- -n_star / 2 * log(2 * pi) - all_det

  eval_one <- function(x, log = default_log) {
    .validate_same_dims(x, dims, "x", reference = "mu")
    x <- array(x, dim = dims)

    xm <- x - mu
    xm_tmp <- xm

    for (k in length(sigmas):1) {
      xm_tmp <- n_prod(xm_tmp, inv_sigmas[[k]], k)
    }

    delta <- sum(as.numeric(xm_tmp) * matrix(c(xm), nrow = n_star))
    loglik <- const - 0.5 * delta

    if (log) loglik else exp(loglik)
  }

  function(x, log = default_log) {
    .eval_density_input(x, eval_one, log)
  }
}

.make_dtskewt_evaluator <- function(mu, skew, sigmas, nu, log = FALSE) {
  default_log <- log
  dims <- .tensor_dims(mu)
  n_star <- prod(dims)

  .validate_same_dims(skew, dims, "skew", reference = "mu")
  sigmas <- .prepare_sigmas(sigmas, dims)

  mu <- array(mu, dim = dims)
  skew <- array(skew, dim = dims)
  inv_sigmas <- lapply(sigmas, invert_safe)
  all_det <- .density_all_det(sigmas, dims)

  skew_tmp <- skew
  for (k in length(sigmas):1) {
    skew_tmp <- n_prod(skew_tmp, inv_sigmas[[k]], k)
  }

  ve_skew <- matrix(c(skew), nrow = n_star)
  rho <- sum(as.numeric(skew_tmp) * ve_skew)
  const <- log(2) + nu / 2 * log(nu / 2) - lgamma(nu / 2) -
    n_star / 2 * log(2 * pi) - all_det

  eval_one <- function(x, log = default_log) {
    .validate_same_dims(x, dims, "x", reference = "mu")
    x <- array(x, dim = dims)

    xm <- x - mu
    xm_tmp <- xm

    for (k in length(sigmas):1) {
      xm_tmp <- n_prod(xm_tmp, inv_sigmas[[k]], k)
    }

    ve_xm <- matrix(c(xm), nrow = n_star)
    xm_tmp <- as.numeric(xm_tmp)
    delta <- sum(xm_tmp * ve_xm)
    xm_skew <- sum(xm_tmp * ve_skew)
    y <- sqrt(rho * (delta + nu))
    log_bessel <- .log_besselK_asym(y, (nu + n_star) / 2)

    loglik <- const + xm_skew -
      (nu + n_star) / 4 * log((delta + nu) / rho) +
      log_bessel - y

    if (log) loglik else exp(loglik)
  }

  function(x, log = default_log) {
    .eval_density_input(x, eval_one, log)
  }
}

.make_dtinvgauss_evaluator <- function(mu, skew, sigmas, kappa, log = FALSE) {
  default_log <- log
  dims <- .tensor_dims(mu)
  n_star <- prod(dims)

  .validate_same_dims(skew, dims, "skew", reference = "mu")
  sigmas <- .prepare_sigmas(sigmas, dims)

  mu <- array(mu, dim = dims)
  skew <- array(skew, dim = dims)
  inv_sigmas <- lapply(sigmas, invert_safe)
  all_det <- .density_all_det(sigmas, dims)

  skew_tmp <- skew
  for (k in seq_along(sigmas)) {
    skew_tmp <- n_prod(skew_tmp, inv_sigmas[[k]], k)
  }

  ve_skew <- matrix(c(skew), nrow = n_star)
  rho <- sum(as.numeric(skew_tmp) * ve_skew)
  const <- log(2) - (n_star + 1) / 2 * log(2 * pi) -
    all_det + kappa

  eval_one <- function(x, log = default_log) {
    .validate_same_dims(x, dims, "x", reference = "mu")
    x <- array(x, dim = dims)

    xm <- x - mu
    xm_tmp <- xm

    for (k in seq_along(sigmas)) {
      xm_tmp <- n_prod(xm_tmp, inv_sigmas[[k]], k)
    }

    ve_xm <- matrix(c(xm), nrow = n_star)
    xm_tmp <- as.numeric(xm_tmp)
    delta <- sum(xm_tmp * ve_xm)
    xm_skew <- sum(xm_tmp * ve_skew)
    y <- sqrt((rho + kappa^2) * (delta + 1))
    log_bessel <- .log_besselK_asym(y, (1 + n_star) / 2)

    loglik <- const + xm_skew -
      (1 + n_star) / 4 * log((delta + 1) / (rho + kappa^2)) +
      log_bessel - y

    if (log) loglik else exp(loglik)
  }

  function(x, log = default_log) {
    .eval_density_input(x, eval_one, log)
  }
}

.make_dtvargamma_evaluator <- function(mu, skew, sigmas, scale, log = FALSE) {
  default_log <- log
  dims <- .tensor_dims(mu)
  n_star <- prod(dims)

  .validate_same_dims(skew, dims, "skew", reference = "mu")
  sigmas <- .prepare_sigmas(sigmas, dims)

  mu <- array(mu, dim = dims)
  skew <- array(skew, dim = dims)
  inv_sigmas <- lapply(sigmas, invert_safe)
  all_det <- .density_all_det(sigmas, dims)

  skew_tmp <- skew
  for (k in length(sigmas):1) {
    skew_tmp <- n_prod(skew_tmp, inv_sigmas[[k]], k)
  }

  ve_skew <- matrix(c(skew), nrow = n_star)
  rho <- sum(as.numeric(skew_tmp) * ve_skew)
  const <- log(2) + scale * log(scale) - n_star / 2 * log(2 * pi) -
    all_det - log(gamma(scale))

  eval_one <- function(x, log = default_log) {
    .validate_same_dims(x, dims, "x", reference = "mu")
    x <- array(x, dim = dims)

    xm <- x - mu
    xm_tmp <- xm

    for (k in length(sigmas):1) {
      xm_tmp <- n_prod(xm_tmp, inv_sigmas[[k]], k)
    }

    ve_xm <- matrix(c(xm), nrow = n_star)
    xm_tmp <- as.numeric(xm_tmp)
    delta <- sum(xm_tmp * ve_xm)
    xm_skew <- sum(xm_tmp * ve_skew)
    y <- sqrt((rho + 2 * scale) * delta)
    log_bessel <- .log_besselK_asym(y, scale - n_star / 2)

    loglik <- const + xm_skew +
      (scale - n_star / 2) / 2 * log(delta / (rho + 2 * scale)) +
      log_bessel - y

    if (log) loglik else exp(loglik)
  }

  function(x, log = default_log) {
    .eval_density_input(x, eval_one, log)
  }
}

.make_dtgenhyper_evaluator <- function(mu, skew, sigmas, lambda, omega, log = FALSE) {
  default_log <- log
  dims <- .tensor_dims(mu)
  n_star <- prod(dims)

  .validate_same_dims(skew, dims, "skew", reference = "mu")
  sigmas <- .prepare_sigmas(sigmas, dims)

  mu <- array(mu, dim = dims)
  skew <- array(skew, dim = dims)
  inv_sigmas <- lapply(sigmas, invert_safe)
  all_det <- .density_all_det(sigmas, dims)

  skew_tmp <- skew
  for (k in seq_along(sigmas)) {
    skew_tmp <- n_prod(skew_tmp, inv_sigmas[[k]], k)
  }

  ve_skew <- matrix(c(skew), nrow = n_star)
  rho <- sum(as.numeric(skew_tmp) * ve_skew)
  Kw <- .log_besselK_asym(omega, lambda) - omega
  const <- -n_star / 2 * log(2 * pi) - all_det - Kw

  eval_one <- function(x, log = default_log) {
    .validate_same_dims(x, dims, "x", reference = "mu")
    x <- array(x, dim = dims)

    xm <- x - mu
    xm_tmp <- xm

    for (k in seq_along(sigmas)) {
      xm_tmp <- n_prod(xm_tmp, inv_sigmas[[k]], k)
    }

    ve_xm <- matrix(c(xm), nrow = n_star)
    xm_tmp <- as.numeric(xm_tmp)
    delta <- sum(xm_tmp * ve_xm)
    xm_skew <- sum(xm_tmp * ve_skew)
    y <- sqrt((rho + omega) * (delta + omega))
    log_bessel <- .log_besselK_asym(y, lambda - n_star / 2)

    loglik <- const + xm_skew +
      (lambda - n_star / 2) / 2 * log((delta + omega) / (rho + omega)) +
      log_bessel - y

    if (log) loglik else exp(loglik)
  }

  function(x, log = default_log) {
    .eval_density_input(x, eval_one, log)
  }
}
