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
rtskewt <- function(n, mu = 0, sigmas = list(matrix(1)), skew = 1, nu = 4) {
  dims <- .tensor_dims(mu)
  .validate_same_dims(skew, dims, "skew", "mu")
  sigmas <- .prepare_sigmas(sigmas, dims)

  # draw tensor variate normals
  tensor_norms <- rtnorm(n, mu = array(0, dim = dims), sigmas)

  # generate inv gamma
  inv_gammas <- rinvgamma(n = n, shape = nu/2, rate = nu/2)

  skewt_draws <- vector("list", n)

  for(i in 1:n) { # scale normals and skew by inv gamma
    skewt_draws[[i]] <- mu + tensor_norms[[i]] * sqrt(inv_gammas[i]) +
                        skew * inv_gammas[i]
  }

  skewt_draws
}
