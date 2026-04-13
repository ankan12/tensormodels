#' dtinvgauss
#'
#' Computes the observed log likelihood for the tensor variate
#' inverse Gaussian.
#'
#' @param x An array of quantiles.
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
  dims <- dim(x)
  o <- length(dims)
  n_star <- prod(dims)

  .validate_same_dims(mu, dims, "mu")
  .validate_same_dims(skew, dims, "skew")
  sigmas <- .prepare_sigmas(sigmas, dims)

  x <- array(x, dim = dims)
  mu <- array(mu, dim = dims)
  skew <- array(skew, dim = dims)

  all_det <- 0

  xm_tmp <- x - mu
  skew_tmp <- skew

  ve_xm <- matrix(c(x - mu), nrow = n_star)
  ve_skew <- matrix(c(skew), nrow = n_star)

  for(k in 1:o) {
    sigd <- sigmas[[k]]

    curr_inv <- invert_safe(sigd)

    xm_tmp <- n_prod(xm_tmp, curr_inv, k)
    skew_tmp <- n_prod(skew_tmp, curr_inv, k)

    all_det <- all_det + (n_star/(2 * nrow(sigd))) * log(det(sigd))
  }

  xm_tmp <- as.numeric(xm_tmp)
  skew_tmp <- as.numeric(skew_tmp)

  rho <- sum(skew_tmp * ve_skew)
  delta <- sum(xm_tmp * ve_xm)
  xm_skew <- sum(xm_tmp * ve_skew)

  y <- sqrt((rho + kappa^2) * (delta + 1))

  loglik <-
        log(2) - (n_star + 1)/2 * log(2 * pi) - all_det + xm_skew +
        kappa - (1 + n_star)/4 * log((delta + 1)/(rho + kappa^2)) +
        log(besselK(y, nu = -(1 + n_star)/2, expon.scaled = TRUE)) - y

  if(log == FALSE) loglik <- exp(loglik)

  loglik
}
