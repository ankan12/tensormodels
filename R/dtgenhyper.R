#' dtgenhyper
#' Computes the density function for the tensor variate generalized hyperbolic.
#'
#' @param x An array of quantiles.
#' @param mu An array of the mean.
#' @param skew An array of the skew.
#' @param sigmas A list of the sigma covariance matrices.
#' @param lambda A parameter that describes the shape of the
#'               generalized inverse Gaussian.
#' @param omega A parameter that describes the shape and scale of the
#'              generalized inverse Gaussian.
#' @param log Defaults to FALSE. If TRUE, return the log of the density.
#'
#' @return The log likelihood.
#'
#' @examples
#' one_draw <- rtgenhyper(n = 1, mu = 2)
#' dtgenhyper(one_draw[[1]], mu = 0, skew = 1,
#'            sigmas = list(matrix(1)), lambda = 2, omega = 2)
#' @export
dtgenhyper <- function(x, mu, skew, sigmas, lambda, omega, log = FALSE) {
  d <- dim(x)
  num_dim <- length(d)
  n_star <- prod(d)

  all_det <- 0
  kroneck_sigmas <- 1

  for(d in length(sigmas):1) {
    kroneck_sigmas <- kronecker(kroneck_sigmas, invert_safe(sigmas[[d]]))
    all_det <- all_det + (n_star/(2 * nrow(sigmas[[d]]))) * log(det(sigmas[[d]]))
  }

  ve_skew <- matrix(c(skew), nrow = n_star)

  loglik <- 0

  centered <- x - mu

  ve_xm <- matrix(c(centered), nrow = n_star)

  xm_skew <- t(ve_xm) %*% kroneck_sigmas %*% ve_skew
  rho <- t(ve_skew) %*% kroneck_sigmas %*% ve_skew
  delta <- t(ve_xm) %*% kroneck_sigmas %*% ve_xm

  Kw <- log(besselK(omega, nu = lambda, expon.scaled = TRUE)) - omega

  loglik <- loglik - (n_star)/2 * log(2 * pi) - all_det - Kw + xm_skew +
    (lambda - n_star/2)/2 * log((delta + omega)/(rho + omega)) +
    log(besselK(sqrt((rho + omega) * (delta + omega)), nu = lambda - n_star/2,
                expon.scaled = TRUE)) - sqrt((rho + omega) * (delta + omega))

  if(log == FALSE) loglik <- exp(loglik)

  loglik
}
