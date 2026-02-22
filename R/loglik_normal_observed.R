#' dtnorm
#'
#' Computes the observed log likelihood for normal.
#'
#' @param draws An array of draws.
#' @param mu An array of the mean.
#' @param sigmas Describes shape of the Generalized Inverse Gaussian.
#'
#' @return The log likelihood.
#'
#' @noRd
loglik_normal_observed <- function(draws, mu, sigmas) {

  n <- length(draws)
  dims <- dim(draws[[1]])
  num_dim <- length(dims)
  n_star <- prod(dims)

  # mu_array <- replicate(n, mu, simplify = "array") |>
  #   aperm(c(num_dim + 1, (1:(num_dim)))) |>
  #   as_tensor()
  #
  # centered <- draws - mu_array

  all_det <- 0
  kroneck_sigmas <- 1

  # for(d in length(sigmas):1) {
  #   kroneck_sigmas <- kronecker(kroneck_sigmas, chol2inv(chol(sigmas[[d]])))
  #   all_det <- all_det + (n_star/(2 * nrow(sigmas[[d]]))) * log(det(sigmas[[d]]))
  # }

  for (d in length(sigmas):1) {
    s_d <- sigmas[[d]]
    U  <- chol(s_d)

    all_det <- all_det + (n_star / (2 * nrow(s_d))) * (2 * sum(log(diag(U))))
    kroneck_sigmas <- kronecker(kroneck_sigmas, chol2inv(U))
  }

  loglik <- 0

  for(i in 1:n) {
    centered <- draws[[i]] - mu

    ve_xm <- matrix(c(centered), nrow = n_star)

    delta <- t(ve_xm) %*% kroneck_sigmas %*% ve_xm

    loglik <- loglik - (n_star)/2 * log(2 * pi) - all_det + (-1/2 * delta)
  }
  loglik
}
