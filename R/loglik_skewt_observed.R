#' loglik_skewt_observed
#'
#' Computes the observed log likelihood for rtskewt for convergence checks.
#'
#' @param draws An array of draws.
#' @param mu An array of the mean.
#' @param skew An array of the skew.
#' @param sigmas A list of the sigma covariance matrices.
#' @param nu Describes df of the invgamma.
#'
#' @return The log likelihood.
#'
#' @noRd
loglik_skewt_observed <- function(draws, mu, skew, sigmas, nu) {
  dims   <- dim(draws)[-1]
  num_dim <- length(dims)
  n      <- dim(draws)[1]
  n_star <- prod(dims)

  mu_array <- replicate(n, mu, simplify = "array") |>
    aperm(c(num_dim + 1, (1:(num_dim))))

  centered <- draws - mu_array

  all_det <- 0
  kroneck_sigmas <- 1

  for(d in length(sigmas):1) {
    kroneck_sigmas <- kronecker(kroneck_sigmas, chol2inv(chol(sigmas[[d]])))
    all_det <- all_det + (n_star/(2 * nrow(sigmas[[d]]))) * log(det(sigmas[[d]]))
  }

  ve_skew <- matrix(c(skew), nrow = n_star)

  loglik <- 0

  for(i in 1:n) {
    ve_xm <- matrix(c(centered[i, , ,]), nrow = n_star)

    xm_skew <- t(ve_xm) %*% kroneck_sigmas %*% ve_skew
    rho <- t(ve_skew) %*% kroneck_sigmas %*% ve_skew
    delta <- t(ve_xm) %*% kroneck_sigmas %*% ve_xm


    loglik <- loglik + log(2) + nu/2 * log(nu/2) -
      log(gamma(nu/2)) - n_star/2 * log(2 * pi) -
      all_det + xm_skew - (nu + n_star)/4 * log((delta + nu)/rho) +
      log(besselK(sqrt(rho * (delta + nu)), nu = -(nu + n_star)/2,
                  expon.scaled = TRUE)) - sqrt(rho * (delta + nu))
  }
  loglik
}
