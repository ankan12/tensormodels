#' tensor_mle_genhyper
#'
#' Computes the EM estimation algorithm for the genhyper distribution.
#' @return A list containing the estimated mean array, skew array, list of
#'   covariance matrices, lambda, and omega.
#'
#' @noRd
tensor_mle_genhyper <- function(data, max_iter = 1000, tol = 1e-6,
                                quiet = TRUE, restrict = NULL) {

  # get dim of input
  n <- length(data)
  dims <- dim(data[[1]])
  num_dim <- length(dims)
  n_star <- prod(dims)

  # dims <- dim(draws)[-1]
  # num_dim <- length(dims)
  # n <- dim(draws)[1]
  # n_star <- prod(dims)

  # Step 1: Initialize vals
  #mu <- apply(X = draws, MARGIN = 2:(num_dim + 1), FUN = mean)
  mu <- simplify2array(data) |> apply(1:num_dim, mean)

  skew <- array(rnorm(prod(dims)), dim = dims)

  #intialize sigmas as identity matrices
  est_sigmas <- lapply(dims, diag)

  logliks <- rep(0, max_iter)

  # different params based on model
  lambda <- 10
  omega <- 10

  for (t in 1:max_iter) {
    # Step 2: Update a, b, c depending on expected values
    skew_compute <- skew
    inv_sigma <- lapply(est_sigmas, invert_safe)

    for (d in seq_along(est_sigmas)) {
      skew_compute <- n_prod(skew_compute, inv_sigma[[d]], d)
    }

    rho <- sum(skew * skew_compute)

    delta_vals <- rep(0, n)

    # mu_array <- replicate(n, mu, simplify = "array") |>
    #   aperm(c(num_dim + 1, (1:(num_dim))))

    #centered <- draws - mu_array

    for (i in 1:n) {
      center_draw <- data[[i]] - mu

      centered_compute <- center_draw

      for (d in seq_along(est_sigmas)) {
        centered_compute <-
          n_prod(centered_compute, inv_sigma[[d]], d)
      }

      delta_vals[i] <- sum(center_draw * centered_compute)
    }

    rho <- rho + omega
    delta_vals <- delta_vals + omega
    param_vals <- lambda - n_star / 2

    k_lambda_1 <- besselK(
      x = sqrt(rho * delta_vals),
      nu = param_vals + 1,
      expon.scaled = TRUE
    )

    k_lambda <- besselK(
      x = sqrt(rho * delta_vals),
      nu = param_vals,
      expon.scaled = TRUE
    )

    a <- sqrt(delta_vals / rho) * (k_lambda_1 / k_lambda)

    b <- sqrt(rho / delta_vals) *
      (k_lambda_1 / k_lambda) -
      (2 * param_vals) / delta_vals

    eps <- 1e-5

    K_plus <- besselK(
      x = sqrt(rho * delta_vals),
      nu = param_vals + eps,
      expon.scaled = TRUE
    )

    K_minus <- besselK(
      x = sqrt(rho * delta_vals),
      nu = param_vals - eps,
      expon.scaled = TRUE
    )

    c <- 1/2 * log(delta_vals / rho) + (log(K_plus) - log(K_minus)) / (2*eps)

    # replace NaN or Inf values
    b[!is.finite(b)] <- 0
    c[!is.finite(c)] <- 0

    # Step 3: Update mu, skew, params
    weight_mean <- mean(a) * b - 1
    weight_skew <- mean(b) - b

    num_mean <- 0
    num_skew <- 0

    for (i in 1:n) {
      num_mean <- num_mean + weight_mean[i] * data[[i]]
      num_skew <- num_skew + weight_skew[i] * data[[i]]
    }

    den_mean <- sum(mean(a) * b) - n

    new_mu <- num_mean / den_mean

    den_skew <- sum(a * mean(b)) - n
    new_skew <- num_skew / den_skew

    new_sigmas <- est_sigmas

    # update params based on model
    update_nu <- function(nu, b, c, n) {
      log(nu / 2) + 1 - digamma(nu / 2) - 1 / n * sum(b + c)
    }
    update_gamma <- function(gamma, a, c) {
      log(gamma) + 1 - digamma(gamma) + mean(c) - mean(a)
    }

    K_plus <- besselK(x = omega, nu = lambda + eps, expon.scaled = TRUE)

    K_minus <- besselK(x = omega, nu = lambda - eps, expon.scaled = TRUE)

    new_lambda <- mean(c) *
      lambda *
      1 /
      ((log(K_plus) - log(K_minus)) / (2 * eps))

    R_lambda <-
      besselK(x = omega, nu = new_lambda + 1, expon.scaled = TRUE) /
      besselK(x = omega, nu = new_lambda, expon.scaled = TRUE)

    R_neg_lambda <-
      besselK(x = omega, nu = -new_lambda + 1, expon.scaled = TRUE) /
      besselK(x = omega, nu = -new_lambda, expon.scaled = TRUE)

    first_deriv <- 1 / 2 * (R_lambda + R_neg_lambda - (mean(a) + mean(b)))

    second_deriv <- 1 /
      2 *
      (R_lambda^2 -
         (1 + 2 * new_lambda) / omega * R_lambda -
         1 +
         R_neg_lambda^2 -
         (1 - 2 * new_lambda) / omega * R_neg_lambda -
         1)

    new_omega <- omega - first_deriv / second_deriv

    scale_prod <- 1

    for (j in 1:num_dim) {
      first <- 0
      second <- 0
      third <- 0
      fourth <- 0

      n_d <- dims[j]

      # covar estimate for mode j
      sigma_j <- matrix(0, nrow = n_d, ncol = n_d)

      # other modes
      other_modes <- (1:num_dim)[-j]

      # Build the whitening operator
      inv_others <- diag(1)

      for (d in rev(other_modes)) {
        inv_others <- kronecker(inv_others, invert_safe(new_sigmas[[d]]))
      }

      A_i <- matricization(new_skew, j)

      for (i in 1:n) {
        Xi_centered <- matricization(data[[i]] - new_mu, j)

        first <- first + b[i] * (Xi_centered %*% inv_others %*% t(Xi_centered))

        second <- second + A_i %*% inv_others %*% t(Xi_centered)

        third <- third + Xi_centered %*% inv_others %*% t(A_i)

        fourth <- fourth + a[i] * A_i %*% inv_others %*% t(A_i)
      }

      sigma_j <- n_d / (n * n_star) * (first - second - third + fourth)

      if (j < num_dim) {
        sigma_j <- sigma_j / (sum(diag(sigma_j))) * n_d
      }

      new_sigmas[[j]] <- sigma_j
    }

    # update all parameters
    mu <- new_mu
    skew <- new_skew
    est_sigmas <- new_sigmas

    lambda <- new_lambda
    omega <- new_omega

    # Step 5: Check convergence

    total_loglik <- 0

    for(i in 1:n) {
      total_loglik <- total_loglik +
        dtgenhyper(data[[i]], mu, skew, est_sigmas, lambda, omega, log = TRUE)
    }

    logliks[t] <- total_loglik

    if(t >= 3) {

      lt_after <- logliks[t]
      lt <- logliks[t - 1]
      lt_before <- logliks[t - 2]

      aitken <- (lt_after - lt)/(lt - lt_before)

      linf <- lt + 1/(1 - aitken) * (lt_after - lt)

      converge <- abs(linf - lt_after)

      if(converge < tol) {
        if(!quiet) message("Converged at iteration ", t)
        break
      }

      if (t %% 50 == 0 & !quiet) {
        cat(sprintf(
          "Iteration %d: mean relative change = %.3e\n",
          t,
          converge
        ))
      }
    }
  }

  if(t == max_iter) message("Reached max iter ", max_iter)

  list(mu = mu, skew = skew, sigmas = est_sigmas,
       lambda = lambda, omega = omega,
       Ew = a, Einvw = b, Elogw = c)
}
