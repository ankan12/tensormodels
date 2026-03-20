#' dtskewt
#'
#' Computes the density function for the tensor variate skewed-t.
#'
#' @param x An array of quantiles.
#' @param mu An array of the mean.
#' @param skew An array of the skew.
#' @param sigmas A list of the sigma covariance matrices.
#' @param nu A parameter that describes the shape and rate of the
#'           inverse gamma.
#' @param log Defaults to FALSE. If TRUE, return the log of the density.
#'
#' @return The log likelihood.
#'
#' @examples
#' one_draw <- rtskewt(n = 1, mu = 2)
#' dtskewt(one_draw[[1]], mu = 0, skew = 1, sigmas = list(matrix(1)), nu = 2)
#' @export
dtskewt <- function(x, mu, skew, sigmas, nu, log = FALSE) {
  d <- dim(x)
  num_dim <- length(d)
  n_star <- prod(d)

  if(is.null(d)) {
    d <- 1
    num_dim <- 1
    n_star <- length(x)
  }

  centered <- x - mu

  all_det <- 0
  kroneck_sigmas <- 1

  for(d in length(sigmas):1) {
    kroneck_sigmas <- kronecker(kroneck_sigmas, invert_safe(sigmas[[d]]))
    all_det <- all_det + (n_star/(2 * nrow(sigmas[[d]]))) *
      determinant(sigmas[[d]], logarithm = TRUE)$modulus[1]
  }

  ve_skew <- matrix(c(skew), nrow = n_star)

  rho <- t(ve_skew) %*% kroneck_sigmas %*% ve_skew

  ve_xm <- matrix(c(x), nrow = n_star)

  xm_skew <- t(ve_xm) %*% kroneck_sigmas %*% ve_skew
  delta <- t(ve_xm) %*% kroneck_sigmas %*% ve_xm

  y <- sqrt(rho * (delta + nu))

  loglik <-
    log(2) + nu/2 * log(nu/2) - lgamma(nu/2) - n_star/2 * log(2 * pi) -
    all_det + xm_skew - (nu + n_star)/4 * log((delta + nu)/rho) +
    log(besselK(y, nu = -(nu + n_star)/2, expon.scaled = TRUE)) - y

  if(log == FALSE) loglik <- exp(loglik)

  loglik
}
