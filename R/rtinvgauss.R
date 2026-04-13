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
#' mean(univar_invgauss)
#' sd(univar_invgauss)
#' @export
#' @importFrom statmod rinvgauss


rtinvgauss <- function(n, mu = 0, sigmas = list(matrix(1)), skew = 1, kappa = 2) {
  dims <- dim(mu)
  .validate_same_dims(skew, dims, "skew", "mu")
  sigmas <- .prepare_sigmas(sigmas, dims)

  # draw tensor variate normals
  tensor_norms <- rtnorm(n, mu = array(0, dim = dims), sigmas)

  # generate inv Gauss
  inv_gauss <- rinvgauss(n = n, mean = 1/kappa, shape = 1)

  invgauss_draws <- vector("list", n)

  for(i in 1:n) { # scale normals and skew by inv Gauss
    invgauss_draws[[i]] <- mu + tensor_norms[[i]] * sqrt(inv_gauss[i]) +
                           skew * inv_gauss[i]
  }

  invgauss_draws
}
