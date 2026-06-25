#' dtnorm
#'
#' Computes the density function for the tensor variate normal.
#'
#' @param x An array of quantiles or a list of arrays.
#' @param mu An array of the mean.
#' @param sigmas A list of covariance matrices.
#' @param log Defaults to FALSE. If TRUE, return the log of the density.
#'
#' @return The log likelihood.
#'
#' @examples
#' dtnorm(x = array(1), mu = 0, sigmas = list(matrix(1)))
#' dnorm(1)
#' dtnorm(x = array(1:12, dim = c(2, 3)), mu = array(0, dim = c(2, 3)),
#'        sigmas = lapply(c(2, 3), diag))
#' @export
dtnorm <- function(x, mu = NULL, sigmas = NULL, log = FALSE) {
  x_is_list <- is.list(x)

  if (x_is_list) {
    if (length(x) == 0) {
      stop("x must be a non-empty list of tensor draws.")
    }

    dims <- dim(x[[1]])

    if (is.null(dims)) {
      stop("Each element of x must be an array.")
    }
  } else {
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
    det_sigd <- det(sigd)

    if (det_sigd == 0) {
      det_sigd <- 1e-8
    }

    all_det <- all_det + (n_star / (2 * nrow(sigd))) * log(det_sigd)
  }

  const <- -n_star / 2 * log(2 * pi) - all_det

  eval_one <- function(x_curr) { # evaluate density for each tensor in list
    .validate_same_dims(x_curr, dims, "x", reference = "mu")

    x_curr <- array(x_curr, dim = dims)

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

  if (!x_is_list) {
    return(eval_one(x))
  }

  vals <- numeric(length(x))

  for (i in seq_along(x)) {
    vals[i] <- eval_one(x[[i]])
  }

  vals
}
