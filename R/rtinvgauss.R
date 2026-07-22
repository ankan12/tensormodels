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
#' @return A `tensor` object containing n draws of the tensor variate
#'   inverse Gaussian distribution, with observations stored on the first axis.
#'
#' @examples
#' univar_invgauss <- rtinvgauss(n = 1e3, mu = 0)
#' mean(univar_invgauss)
#' sd(univar_invgauss)
#' @export
#' @importFrom statmod rinvgauss


rtinvgauss <- function(n, mu = 0, sigmas = list(matrix(1)), skew = 1, kappa = 2) {
  dims <- .tensor_dims(mu)
  .validate_same_dims(skew, dims, "skew", "mu")
  sigmas <- .prepare_sigmas(sigmas, dims)

  # draw tensor variate normals; observations are stored on the first axis
  tensor_norms <- unclass(rtnorm(n, mu = array(0, dim = dims), sigmas))

  # generate inv Gauss
  inv_gauss <- rinvgauss(n = n, mean = 1/kappa, shape = 1)

  mu_array <- array(rep(as.numeric(mu), each = n), dim = c(n, dims))
  skew_array <- array(rep(as.numeric(skew), each = n), dim = c(n, dims))

  draws <- sweep(tensor_norms, 1L, sqrt(inv_gauss), "*") +
    mu_array +
    sweep(skew_array, 1L, inv_gauss, "*")

  tensor(draws, obs = 1L)
}
