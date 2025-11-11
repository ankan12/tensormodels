#' rtinvgass
#'
#' Simulate random draws from the tensor variate inverse Gaussian distribution.
#'
#' @param n An integer stating the sample size.
#' @param mu An array containing the mean values of and dims of each draw.
#' @param sigmas A list of covariance matrices. Defaults to the identity.
#' @param skew An array that determines how skewed the distribution is.
#' @param kappa Shape for inverse Gaussian distribution.
#'
#' @return An array containing n draws of the tensor variate inverse Gaussian distribution.
#'
#' @examples
#' univar_invgauss <- rtinvgauss(n = 1e3, mu = 0)
#' mean(univar_skewt)
#' sd(univar_skewt)
#' @export
#' @importFrom statmod rinvgauss


rtinvgauss <- function(n, mu = 0, sigmas = 1, skew = 1, kappa = 2) {
  dims <- dim(mu)

  # mu was a scalar
  if(is.vector(mu)) dims <- 1

  # draw tensor variate normals
  tensor_norms <- rtnorm(n, mu = mu, sigmas)

  all_dim <- dim(tensor_norms)

  # generate inv Gauss
  inv_gauss <- rinvgauss(n = n, mean = 1, shape = kappa)

  # scale normals by inv Gauss
  scale_norms <- tensor_norms * array(sqrt(inv_gauss), dim = all_dim)

  # scale skew by inv Gauss
  scale_skew <- array(rep(skew, each = n), dim = all_dim) * array(inv_gauss, dim = all_dim)

  array(rep(mu, each = n), dim = c(n, dims)) + scale_skew + scale_norms
}
