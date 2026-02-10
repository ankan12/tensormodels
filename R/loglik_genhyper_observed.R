#' loglik_genhyper_observed
#'
#' Computes the observed log likelihood for genhyper for convergence checks.
#'
#' @param draws An array of draws.
#' @param mu An array of the mean.
#' @param skew An array of the skew.
#' @param sigmas Describes shape of the Generalized Inverse Gaussian.
#' @param omega Describes shape/scale of the Generalized Inverse Gaussian.
#'
#' @return The log likelihood.
#'
#' @noRd
loglik_genhyper_observed <- function(draws, mu, skew, sigmas, lambda, omega) {
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

    Kw <- log(besselK(omega, nu = lambda, expon.scaled = TRUE)) - omega

    loglik <- loglik - (n_star)/2 * log(2 * pi) - all_det - Kw + xm_skew +
      (lambda - n_star/2)/2 * log((delta + omega)/(rho + omega)) +
      log(besselK(sqrt((rho + omega) * (delta + omega)), nu = lambda - n_star/2,
                  expon.scaled = TRUE)) - sqrt((rho + omega) * (delta + omega))
  }
  loglik
}
