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

rtvargamma <- function(n, mu = 0, sigmas = list(matrix(1)), skew = 1, scale = 2) {
  dims <- dim(mu)
  .validate_same_dims(skew, dims, "skew", "mu")
  sigmas <- .prepare_sigmas(sigmas, dims)

  # draw tensor variate normals
  tensor_norms <- rtnorm(n, mu = array(0, dim = dims), sigmas)

  # generate gamma draws
  gammas <- rgamma(n = n, shape = scale, rate = scale)

  vargamma_draws <- vector("list", n)

  for(i in 1:n) { # scale normals and skew by gamma
    vargamma_draws[[i]] <- mu + tensor_norms[[i]] * sqrt(gammas[i]) + skew * gammas[i]
  }

  vargamma_draws
}
