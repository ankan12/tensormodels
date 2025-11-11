#' rtvargamma
#'
#' Simulate random draws from the tensor variate variance gamma distribution.
#'
#' @param n An integer stating the sample size.
#' @param mu An array containing the mean values of and dims of each draw.
#' @param sigmas A list of covariance matrices. Defaults to the identity.
#' @param skew An array that determines how skewed the distribution is.
#' @param scale Scale parameter for gamma.
#'
#' @return An array containing n draws of the tensor variate variance gamma distribution.
#'
#' @examples
#' univar_vargam <- rtvargamma(n = 1e3, mu = 0)
#' mean(univar_vargam)
#' sd(univar_vargam)
#' @export

rtvargamma <- function(n, mu = 0, sigmas = 1, skew = 1, scale = 2) {
  dims <- dim(mu)

  # mu was a scalar
  if(is.vector(mu)) dims <- 1

  # draw tensor variate normals
  tensor_norms <- rtnorm(n, mu = array(0, dim = dims), sigmas)

  all_dim <- dim(tensor_norms)

  # generate gamma draws
  gammas <- rgamma(n = n, shape = scale, rate = scale)

  # scale normals by gamma
  scale_norms <- tensor_norms * array(sqrt(gammas), dim = all_dim)

  # scale skew by gamma
  scale_skew <- array(rep(skew, each = n), dim = all_dim) * array(gammas, dim = all_dim)

  array(rep(mu, each = n), dim = c(n, dims)) + scale_skew + scale_norms
}
