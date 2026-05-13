#' rtgenhyper
#'
#' Simulate random draws from the tensor variate generalized hyperbolic distribution.
#'
#' @param n An integer stating the sample size.
#' @param mu An array containing the mean values of and dims of each draw.
#' @param skew An array that determines how skewed the distribution is.
#' @param sigmas A list of covariance matrices. Defaults to the identity.
#' @param lambda Shape/scale parameter for generalized inverse Gaussian distribution.
#' @param omega Shape parameter for generalized inverse Gaussian distribution.
#'
#' @return An array containing n draws of the tensor generalized hyperbolic distribution.
#'
#' @examples
#' univar_genhyper <- rtgenhyper(n = 1e3, mu = 0)
#' mean(univar_genhyper)
#' sd(univar_genhyper)
#' @export
#' @importFrom GIGrvg rgig


rtgenhyper <- function(n, mu = 0, skew = 1, sigmas = list(matrix(1)),
                       lambda = 2, omega = 2) {
  dims <- .tensor_dims(mu)
  .validate_same_dims(skew, dims, "skew", "mu")
  sigmas <- .prepare_sigmas(sigmas, dims)

  # draw tensor variate normals
  tensor_norms <- rtnorm(n, mu = array(0, dim = dims), sigmas)

  # generate inv GIG
  inv_gig <- rgig(n, lambda = lambda, chi = omega, psi = omega)

  genhyper_draws <- vector("list", n)

  for(i in 1:n) { # scale normals and skew by inv GIG
    genhyper_draws[[i]] <- mu + tensor_norms[[i]] * sqrt(inv_gig[i]) + skew * inv_gig[i]
  }

  genhyper_draws
}
