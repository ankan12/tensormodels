#' dtinvgauss
#'
#' Computes the observed log likelihood for the tensor variate
#' inverse Gaussian.
#'
#' @param x An array of quantiles or a list of arrays.
#' @param mu An array of the mean.
#' @param skew An array of the skew.
#' @param sigmas A list of the sigma covariance matrices.
#' @param kappa A parameter that describes the shape of the
#'           inverse Gaussian
#' @param log Defaults to FALSE. If TRUE, return the log of the density.
#'
#' @return The log likelihood.
#'
#' @examples
#' dtinvgauss(array(1), mu = 0, skew = array(1), sigmas = list(matrix(1)), kappa = 2)
#' @export
dtinvgauss <- function(x, mu, skew, sigmas, kappa, log = FALSE) {
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
  .validate_same_dims(skew, dims, "skew")
  sigmas <- .prepare_sigmas(sigmas, dims)

  mu <- array(mu, dim = dims)
  skew <- array(skew, dim = dims)

  all_det <- 0
  skew_tmp <- skew

  ve_skew <- matrix(c(skew), nrow = n_star)
  inv_sigmas <- lapply(sigmas, invert_safe)

  for(k in 1:o) {
    sigd <- sigmas[[k]]

    skew_tmp <- n_prod(skew_tmp, inv_sigmas[[k]], k)

    all_det <- all_det + (n_star/(2 * nrow(sigd))) * log(det(sigd))
  }

  skew_tmp <- as.numeric(skew_tmp)

  rho <- sum(skew_tmp * ve_skew)

  eval_one <- function(x_curr) {
    .validate_same_dims(x_curr, dims, "x", reference = "mu")

    x_curr <- array(x_curr, dim = dims)

    xm <- x_curr - mu
    xm_tmp <- xm
    ve_xm <- matrix(c(xm), nrow = n_star)

    for(k in 1:o) {
      xm_tmp <- n_prod(xm_tmp, inv_sigmas[[k]], k)
    }

    xm_tmp <- as.numeric(xm_tmp)

    delta <- sum(xm_tmp * ve_xm)
    xm_skew <- sum(xm_tmp * ve_skew)

    y <- sqrt((rho + kappa^2) * (delta + 1))
    log_bessel <- .log_besselK_asym(y, (1 + n_star) / 2)

    loglik <-
          log(2) - (n_star + 1)/2 * log(2 * pi) - all_det + xm_skew +
          kappa - (1 + n_star)/4 * log((delta + 1)/(rho + kappa^2)) +
          log_bessel - y

    if(log == FALSE) loglik <- exp(loglik)

    loglik
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
