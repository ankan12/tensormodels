#' rtgenhyper
#'
#' Simulate random draws from the tensor variate generalized hyperbolic distribution.
#'
#' @param n An integer stating the sample size.
#' @param mu An array containing the mean values of and dims of each draw.
#' @param sigmas A list of covariance matrices. Defaults to the identity.
#' @param skew An array that determines how skewed the distribution is.
#' @param omega Shape parameter for generalized inverse Gaussian distribution.
#' @param lambda Shape/scale parameter for generalized inverse Gaussian distribution.
#'
#' @return An array containing n draws of the tensor generalized hyperbolic distribution.
#'
#' @examples
#' univar_skewt <- rtskewt(n = 1e3, mu = 0, skew = 0, v = 4)
#' mean(univar_skewt)
#' sd(univar_skewt)
#' @export
#' @importFrom GIGrvg rgig


rtgenhyper <- function(n = 50, mu = 0, sigmas = NULL, skew = 1, omega = 2, lambda = 2) {
  all_dims <- dim(mu)

  # mu was a scalar
  if(is.null(all_dims)) all_dims <- 1

  # draw tensor variate normals
  tensor_norms <- rtnorm(n, mu = mu, sigmas)

  # generate inv GIG
  inv_gig <- rgig(n = n, lambda = lambda, chi = omega, psi = 2)

  # scale normals by inv GIG
  scale_norms <- sweep(tensor_norms, 1, sqrt(inv_gig), `*`)

  # scale skew by inv GIG
  scale_skew <- sweep(array(rep(skew, each = n), dim = c(n, all_dims)), 1, inv_gig, `*`)

  array(rep(mu, each = n), dim = c(n, all_dims)) + scale_skew + scale_norms
}
