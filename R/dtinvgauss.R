#' dtinvgauss
#'
#' Computes the observed log likelihood for the tensor variate
#' inverse Gaussian.
#'
#' @param x An array of quantiles.
#' @param mu An array of the mean.
#' @param skew An array of the skew.
#' @param sigmas A list of the sigma covariance matrices.
#' @param nu A parameter that describes the shape of the
#'           inverse Gaussian
#' @param log Defaults to FALSE. If TRUE, return the log of the density.
#'
#' @return The log likelihood.
#'
#'#' @examples
#' one_draw <- rtinvgauss(n = 1, mu = 2)
#' dtinvgauss(one_draw[[1]], mu = 0, skew = 1,
#'            sigmas = list(matrix(1)), kappa = 2)
#' @export
dtinvgauss <- function(x, mu, skew, sigmas, kappa, log = FALSE) {
  d <- dim(x)
  num_dim <- length(d)
  n_star <- prod(d)

  all_det <- 0
  kroneck_sigmas <- 1

  for(d in length(sigmas):1) {
    kroneck_sigmas <- kronecker(kroneck_sigmas, invert_safe(sigmas[[d]]))
    all_det <- all_det + (n_star/(2 * nrow(sigmas[[d]]))) *
               determinant(sigmas[[d]], logarithm = TRUE)$modulus[1]
  }

  ve_skew <- matrix(c(skew), nrow = n_star)

  rho <- t(ve_skew) %*% kroneck_sigmas %*% ve_skew

  loglik <- 0

  centered <- x - mu

  ve_xm <- matrix(c(centered), nrow = n_star)

  xm_skew <- t(ve_xm) %*% kroneck_sigmas %*% ve_skew
  delta <- t(ve_xm) %*% kroneck_sigmas %*% ve_xm

  y <- sqrt((rho + kappa^2) * (delta + 1))

  loglik <- loglik + log(2) - (n_star + 1)/2 * log(2 * pi) - all_det +
              xm_skew + kappa - (1 + n_star)/4 * log((delta + 1)/(rho + kappa^2)) +
              log(besselK(y, nu = -(1 + n_star)/2, expon.scaled = TRUE)) - y

  if(log == FALSE) loglik <- exp(loglik)

  loglik
}
