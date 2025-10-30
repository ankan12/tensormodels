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
#' univar_skewt <- rtskewt(n = 1e3, mu = 0, skew = 0, v = 4)
#' mean(univar_skewt)
#' sd(univar_skewt)
#' @export
#' @importFrom invgamma rinvgamma


rtskewt <- function(n, mu = 0, sigmas = 1, skew = 1, nu = 2) {
  dims <- dim(mu)

  # mu was a scalar
  if(is.vector(mu)) dims <- 1

  # draw tensor variate normals
  tensor_norms <- rtnorm(n, mu = mu, sigmas)

  # generate inv gamma
  inv_gammas <- rinvgamma(n = n, shape = nu/2, rate = nu/2)

  # scale normals by inv gamma
  scale_norms <- sweep(tensor_norms, 1, sqrt(inv_gammas), `*`)

  # scale skew by inv gamma
  scale_skew <- sweep(array(rep(skew, each = n), dim = c(n, dims)), 1, inv_gammas, `*`)

  array(rep(mu, each = n), dim = c(n, dims)) + scale_skew + scale_norms
}
