#' dtnorm
#'
#' Computes the density function for the tensor variate normal.
#'
#' @param x A numeric vector/array, a `tensor`, or a `tensor` object. A
#'   sample returns one density per observation.
#' @param mu An array of the mean.
#' @param sigmas A list of covariance matrices.
#' @param log Defaults to FALSE. If TRUE, return the log of the density.
#'
#' @return A single density for one tensor, or a numeric vector of densities
#'   for a tensor object.
#'
#' @examples
#' dtnorm(x = array(1), mu = 0, sigmas = list(matrix(1)))
#' dnorm(1)
#' dtnorm(x = array(1:12, dim = c(2, 3)), mu = array(0, dim = c(2, 3)),
#'        sigmas = lapply(c(2, 3), diag))
#' @export
dtnorm <- function(x, mu = NULL, sigmas = NULL, log = FALSE) {
  x_is_tensor <- inherits(x, "tensor")

  if (is.list(x)) {
    stop(
      "`x` must be a tensor/array or a `tensor` object; lists are not supported.",
      call. = FALSE
    )
  }

  if (x_is_tensor) {
    if (n_draws(x) == 0L) {
      stop("x must contain at least one tensor observation.", call. = FALSE)
    }

    dims <- draw_shape(x)
  } else {
    if (!is.numeric(x)) {
      stop(
        "`x` must be numeric, an array, a `tensor`, or a `tensor` object.",
        call. = FALSE
      )
    }

    dims <- .tensor_dims(x)
  }

  o <- length(dims)
  n_star <- prod(dims)

  .validate_same_dims(mu, dims, "mu")
  sigmas <- .prepare_sigmas(sigmas, dims)

  inv_sigmas <- lapply(sigmas, invert_safe)

  mu <- array(mu, dim = dims)
  all_det <- 0

  for (k in seq_along(sigmas)) {
    sigd <- sigmas[[k]]
    logdet_sigd <- .logdet_safe(sigd, eps = 1e-10, warn = FALSE)
    all_det <- all_det + (n_star / (2 * nrow(sigd))) * logdet_sigd
  }

  const <- -n_star / 2 * log(2 * pi) - all_det

  eval_one <- function(x_curr) {
    .validate_same_dims(x_curr, dims, "x", reference = "mu")

    if (inherits(x_curr, "tensor")) {
      x_curr <- .tensor_single_draw_array(x_curr)
    } else {
      x_curr <- array(x_curr, dim = dims)
    }

    xm <- x_curr - mu
    xm_tmp <- xm
    ve_xm <- matrix(c(xm), nrow = n_star)

    for (k in length(sigmas):1) {
      xm_tmp <- n_prod(xm_tmp, inv_sigmas[[k]], k)
    }

    delta <- sum(as.numeric(xm_tmp) * ve_xm)

    loglik <- const - 0.5 * delta

    if (log) {
      loglik
    } else {
      exp(loglik)
    }
  }

  if (!x_is_tensor) {
    return(eval_one(x))
  }

  vals <- numeric(n_draws(x))

  for (i in seq_len(n_draws(x))) {
    vals[i] <- eval_one(pull_draw(x, i))
  }

  vals
}
