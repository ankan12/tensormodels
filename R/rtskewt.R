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


rtskewt <- function(n = 50, mu = 0, sigmas = NULL, skew = 1, v = 2) {
  all_dims <- dim(mu)

  # mu was a scalar
  if(is.null(all_dims)) all_dims <- 1

  # draw tensor variate normals
  tensor_norms <- rtnorm(n, mu = mu, sigmas)

  # generate inv gamma
  inv_gammas <- rinvgamma(n = n, shape = v/2, rate = v/2)

  # scale normals by inv gamma
  scale_norms <- sweep(tensor_norms, 1, sqrt(inv_gammas), `*`)

  # scale skew by inv gamma
  scale_skew <- sweep(array(rep(skew, each = n), dim = c(n, all_dims)), 1, inv_gammas, `*`)

  array(rep(mu, each = n), dim = c(n, all_dims)) + scale_skew + scale_norms
}
