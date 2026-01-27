#' loglik_invgauss_observed
#'
#' Computes the observed log likelihood for invgauss for convergence checks.
#'
#' @param draws An array of draws.
#' @param mu An array of the mean.
#' @param skew An array of the skew.
#' @param sigmas A list of the sigma covariance matrices.
#' @param kappa Describes shape of the inverse Gaussian.
#'
#' @return The log likelihood.
#'
#' @noRd
loglik_invgauss_observed <- function(draws, mu, skew, sigmas, kappa) {
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
    all_det <- all_det + (n_star/(2 * nrow(sigmas[[d]]))) *
               determinant(sigmas[[d]], logarithm = TRUE)$modulus[1]
  }

  ve_skew <- matrix(c(skew), nrow = n_star)

  rho <- t(ve_skew) %*% kroneck_sigmas %*% ve_skew

  loglik <- 0

  for(i in 1:n) {
    ve_xm <- matrix(c(centered[i, , ,]), nrow = n_star)

    xm_skew <- t(ve_xm) %*% kroneck_sigmas %*% ve_skew
    delta <- t(ve_xm) %*% kroneck_sigmas %*% ve_xm

    x <- sqrt((rho + kappa^2) * (delta + 1))

    loglik <- loglik + log(2) - (n_star + 1)/2 * log(2 * pi) - all_det +
                xm_skew + kappa - (1 + n_star)/4 * log((delta + 1)/(rho + kappa^2)) +
                log(besselK(x, nu = -(1 + n_star)/2, expon.scaled = TRUE)) - x
  }
  loglik
}
