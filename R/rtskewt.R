#' rtskewt
#'
#' Simulate random draws from the tensor variate skewed t distribution.
#'
#' @param n An integer stating the sample size.
#' @param mu An array containing the mean values of and dims of each draw.
#' @param sigmas A list of covariance matrices. Defaults to the identity.
#' @param skew An array that determines how skewed the distribution is.
#' @param v Degrees of freedom
#'
#' @return An array containing n draws of the tensor variate skewed t distribution.
#'
#' @examples
#' univar_skewt <- rtskewt(n = 1e3, mu = 0)
#' mean(univar_skewt)
#' sd(univar_skewt)
#' @export
#' @importFrom invgamma rinvgamma
rtskewt <- function(n, mu = 0, sigmas = 1, skew = 1, nu = 4) {
  dims <- dim(mu)

  # mu was a scalar
  if(is.vector(mu)) dims <- 1

  # draw tensor variate normals
  tensor_norms <- rtnorm(n, mu = array(0, dim = dims), sigmas)

  all_dim <- dim(tensor_norms)

  # generate inv gamma
  inv_gammas <- rinvgamma(n = n, shape = nu/2, rate = nu/2)

  # scale normals by inv gamma
  scale_norms <- tensor_norms * array(sqrt(inv_gammas), dim = all_dim)

  # scale skew by inv gamma
  scale_skew <- array(rep(skew, each = n), dim = all_dim) * array(inv_gammas, dim = all_dim)

  array(rep(mu, each = n), dim = all_dim) + scale_skew + scale_norms
}
