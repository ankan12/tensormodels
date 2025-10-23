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
#' univar_vargam <- rtskewt(n = 1e3, mu = 0, skew = 0, scale = 2)
#' mean(univar_vargam)
#' sd(univar_vargam)
#' @export

rtvargamma <- function(n = 50, mu = 0, sigmas = NULL, skew = 1, scale = 2) {
  all_dims <- dim(mu)

  # mu was a scalar
  if(is.null(all_dims)) all_dims <- 1

  # draw tensor variate normals
  tensor_norms <- rtnorm(n, mu = mu, sigmas)

  gammas <- rgamma(n = n, shape = scale, rate = scale)

  # scale normals by gamma
  scale_norms <- sweep(tensor_norms, 1, sqrt(gammas), `*`)

  # scale skew by gamma
  scale_skew <- sweep(array(rep(skew, each = n), dim = c(n, all_dims)), 1, gammas, `*`)

  array(rep(mu, each = n), dim = c(n, all_dims)) + scale_skew + scale_norms
}
