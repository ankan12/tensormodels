#' dtvargamma
#'
#' Computes the density function for the tensor variate gamma.
#'
#' @param x An array of quantiles.
#' @param mu An array of the mean.
#' @param skew An array of the skew.
#' @param sigmas A list of the sigma covariance matrices.
#' @param scale A parameter that describes the shape and rate of the gamma.
#' @param log Defaults to FALSE. If TRUE, return the log of the density.
#'
#' @return The log likelihood.
#'
#' @examples
#' one_draw <- rtvargamma(n = 1, mu = 2)
#' dtskewt(one_draw[[1]], mu = 0, skew = 1, sigmas = list(matrix(1)), scale = 4)
#' @export
dtvargamma <- function(x, mu, skew, sigmas, scale, log = FALSE) {
  d <- dim(x)
  num_dim <- length(d)
  n_star <- prod(d)

  all_det <- 0
  kroneck_sigmas <- 1

  for(d in length(sigmas):1) {
    kroneck_sigmas <- kronecker(kroneck_sigmas, chol2inv(chol(sigmas[[d]])))
    all_det <- all_det + (n_star/(2 * nrow(sigmas[[d]]))) * log(det(sigmas[[d]]))
  }

  ve_skew <- matrix(c(skew), nrow = n_star)

  loglik <- 0

  centered <- x - mu

  ve_xm <- matrix(c(centered), nrow = n_star)

  xm_skew <- t(ve_xm) %*% kroneck_sigmas %*% ve_skew
  rho <- t(ve_skew) %*% kroneck_sigmas %*% ve_skew
  delta <- t(ve_xm) %*% kroneck_sigmas %*% ve_xm

  loglik <- loglik + log(2) + scale * log(scale) - (n_star)/2 * log(2 * pi) -
    all_det - log(gamma(scale)) + xm_skew +
    (scale - (n_star/2))/2 * log(delta/(rho + 2 * scale))
    + log(besselK(sqrt((rho + 2 * scale) * (delta)), nu = scale - n_star/2,
                expon.scaled = TRUE)) - sqrt((rho + 2 * scale) * (delta))

  if(log == FALSE) loglik <- exp(loglik)

  loglik
}
