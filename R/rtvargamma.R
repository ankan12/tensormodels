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
#' @return A `tensor` object containing n draws of the tensor variate
#'   variance gamma distribution, with observations stored on the first axis.
#'
#' @examples
#' univar_vargam <- rtvargamma(n = 1e3, mu = 0)
#' mean(univar_vargam)
#' sd(univar_vargam)
#' @export

rtvargamma <- function(n, mu = 0, sigmas = list(matrix(1)), skew = 1, scale = 2) {
  dims <- .tensor_dims(mu)
  .validate_same_dims(skew, dims, "skew", "mu")
  sigmas <- .prepare_sigmas(sigmas, dims)

  # draw tensor variate normals; observations are stored on the first axis
  tensor_norms <- unclass(rtnorm(n, mu = array(0, dim = dims), sigmas))

  # generate gamma draws
  gammas <- rgamma(n = n, shape = scale, rate = scale)

  mu_array <- array(rep(as.numeric(mu), each = n), dim = c(n, dims))
  skew_array <- array(rep(as.numeric(skew), each = n), dim = c(n, dims))

  draws <- sweep(tensor_norms, 1L, sqrt(gammas), "*") +
    mu_array +
    sweep(skew_array, 1L, gammas, "*")

  tensor(draws, obs = 1L)
}
