#' tensor_mle_invgauss
#'
#' Computes the EM estimation algorithm for the invgamma distribution.
#' @return A list containing the estimated mean array, skew array, list of
#'   covariance matrices, and kappa,
#'
#' @noRd
tensor_mle_invgauss <- function(data, max_iter, tol,
                                quiet = TRUE, restrict = NULL) {

  n <- length(data)
  dims <- dim(data[[1]])
  num_dim <- length(dims)
  n_star <- prod(dims)

  # dims <- dim(draws)[-1]
  # num_dim <- length(dims)
  # n <- dim(draws)[1]
  # n_star <- prod(dims)

  # Step 1: Initialize vals
  #mu <- apply(X = data, MARGIN = 2:(num_dim + 1), FUN = mean)
  mu <- simplify2array(data) |> apply(1:num_dim, mean)

  skew <- array(rnorm(prod(dims)), dim = dims)

  logliks <- rep(0, max_iter)

  est_sigmas <- lapply(dims, diag)

  for(k in 1:num_dim) {
    tot_sum <- 0
    for(i in 1:n) {
      curr_unfold <- matricization(data[[i]] - mu, k)

      tot_sum <- tot_sum + (curr_unfold %*% t(curr_unfold))
    }

    tot_sum <- tot_sum * dims[k]/(n * n_star)

    tot_sum <- tot_sum/(sum(diag(tot_sum))) * dims[k]

    est_sigmas[[k]] <- tot_sum
  }

  # different params based on model
  kappa <- 2

  for (t in 1:max_iter) {
    # Step 2: Update a, b, c depending on expected values
    skew_compute <- skew
    inv_sigma <- lapply(est_sigmas, invert_safe)

    for (d in seq_along(est_sigmas)) {
      skew_compute <- n_prod(skew_compute, inv_sigma[[d]], d)
    }

    rho <- sum(skew * skew_compute)

    delta_vals <- rep(0, n)

    for (i in 1:n) {
      center_draw <- data[[i]] - mu

      centered_compute <- center_draw

      for (d in seq_along(est_sigmas)) {
        centered_compute <- n_prod(centered_compute, inv_sigma[[d]], d)
      }
      delta_vals[i] <- sum(center_draw * centered_compute)
    }

    rho <- rho + kappa^2
    delta_vals <- delta_vals + 1
    param_vals <- -(1 + n_star) / 2

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

    new_kappa = n / (sum(a))

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

      scale_prod <- 1

      if (j < num_dim) {
        c_d <- sigma_j[1, 1]

        if (!is.finite(c_d) || c_d == 0) c_d <- 1

        scale_prod <- scale_prod * c_d

        sigma_j <- sigma_j / c_d
      }

      new_sigmas[[j]] <- sigma_j
    }

    new_sigmas[[num_dim]] <- new_sigmas[[num_dim]] * scale_prod

    # update all parameters
    mu <- new_mu
    skew <- new_skew
    est_sigmas <- new_sigmas
    kappa <- new_kappa

    # Step 5: Check convergence

    total_loglik <- 0

    for(i in 1:n) {
      total_loglik <- total_loglik +
        dtinvgauss(data[[i]], mu, skew, est_sigmas, kappa, log = TRUE)
    }

    logliks[t] <- total_loglik

    if(t >= 3) {

      lt_after <- logliks[t]
      lt <- logliks[t - 1]
      lt_before <- logliks[t - 2]

      converge <- logliks[t] - logliks[t-1]

      if(converge < tol) {
        if(!quiet) message("Converged at iteration ", t)
        break
      }

      if (t %% 50 == 0 & !quiet) {
        cat(sprintf(
          "Iteration %d: criterion based on Aitken = %.3e\n",
          t,
          converge
        ))
      }
    }
  }

  k <- n_star + sum((dims * (dims+1))/2) - (num_dim - 1) + 1

  list(mu = mu, skew = skew, sigmas = est_sigmas, kappa = kappa,
       Ew = a, Einvw = b, Elogw = c,
       loglik = logliks[t], BIC = k * log(n) - 2 * logliks[t])
}
