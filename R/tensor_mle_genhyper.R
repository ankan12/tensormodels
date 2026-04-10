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
  o <- length(dims)
  n_star <- prod(dims)

  # Step 1: Initialize vals

  # different params based on model
  lambda <- 2
  omega <- 2

  R_lambda <- besselK(x = omega, nu = lambda + 1, expon.scaled = TRUE) /
              besselK(x = omega, nu = lambda, expon.scaled = TRUE)

  flat_draws <- simplify2array(data)

  # E[X] = M + K_{lambda+1}(omega)/K_{lambda}(omega) * skew
  mean_draws <- apply(flat_draws, 1:o, mean)
  median_draws <- apply(flat_draws, 1:o, median)

  mu <- median_draws

  skew <- (mean_draws - median_draws)/R_lambda

  sigmas <- lapply(dims, diag)

  logliks <- rep(0, max_iter)

  for(k in 1:o) {
    tot_sum <- 0
    for(i in 1:n) {
      curr_unfold <- matricization(data[[i]] - mu, k)

      tot_sum <- tot_sum + tcrossprod(curr_unfold)
    }

    tot_sum <- tot_sum * dims[k]/(n * n_star)

    tot_sum <- tot_sum/(sum(diag(tot_sum))) * dims[k]

    sigmas[[k]] <- tot_sum
  }

  for (t in 1:max_iter) {
    # Step 2: Update a, b, c depending on expected values
    skew_compute <- skew
    inv_sigma <- lapply(sigmas, invert_safe)

    for (d in seq_along(sigmas)) {
      skew_compute <- n_prod(skew_compute, inv_sigma[[d]], d)
    }

    rho <- sum(skew * skew_compute)

    delta_vals <- rep(0, n)

    for (i in 1:n) {
      center_draw <- data[[i]] - mu

      centered_compute <- center_draw

      for (d in seq_along(sigmas)) {
        centered_compute <- n_prod(centered_compute, inv_sigma[[d]], d)
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

    new_sigmas <- sigmas

    # update params based on model
    update_nu <- function(nu, b, c, n) {
      log(nu / 2) + 1 - digamma(nu / 2) - 1 / n * sum(b + c)
    }

    update_gamma <- function(gamma, a, c) {
      log(gamma) + 1 - digamma(gamma) + mean(c) - mean(a)
    }

    K_plus <- besselK(x = omega, nu = lambda + eps, expon.scaled = TRUE)

    K_minus <- besselK(x = omega, nu = lambda - eps, expon.scaled = TRUE)

    new_lambda <- mean(c) * lambda *
                  1 / ((log(K_plus) - log(K_minus)) / (2 * eps))

    R_lambda <-
      besselK(x = omega, nu = new_lambda + 1, expon.scaled = TRUE) /
      besselK(x = omega, nu = new_lambda, expon.scaled = TRUE)

    R_neg_lambda <-
      besselK(x = omega, nu = -new_lambda + 1, expon.scaled = TRUE) /
      besselK(x = omega, nu = -new_lambda, expon.scaled = TRUE)

    first_deriv <- 1 / 2 * (R_lambda + R_neg_lambda - (mean(a) + mean(b)))

    second_deriv <- 1 /2 * (R_lambda^2 - (1 + 2 * new_lambda) / omega * R_lambda -
                            1 + R_neg_lambda^2 -
                            (1 - 2 * new_lambda) / omega * R_neg_lambda - 1)

    new_omega <- omega - first_deriv / second_deriv

    for (j in 1:o) {
      inv_new_sigma <- lapply(new_sigmas, invert_safe) # compute sigmas

      n_d <- dims[j]
      other_modes <- (1:o)[-j]

      first <- matrix(0, n_d, n_d)
      x_sum  <- matrix(0, n_d, prod(dims[-j]))

      skew_tmp <- new_skew # compute skew once for each dim

      for (d in other_modes) { # multiply by inverse covars
        skew_tmp <- n_prod(skew_tmp, inv_new_sigma[[d]], d)
      }

      flat_skew_tmp <- matricization(skew_tmp, j) # flatten the skews
      flat_skew <- matricization(new_skew, j)

      for (i in 1:n) {
        xm <- data[[i]] - new_mu # take centered draw
        flat_xm <- matricization(xm, j)

        x_sum <- x_sum + flat_xm

        xm_tmp <- xm

        for (d in other_modes) { # multiply by inverse covars
          xm_tmp <- n_prod(xm_tmp, inv_new_sigma[[d]], d)
        }

        flat_xm_tmp <- matricization(xm_tmp, j)

        first <- first + b[i] * (flat_xm_tmp %*% t(flat_xm))
      }

      second <- flat_skew_tmp %*% t(x_sum)
      third  <- x_sum %*% t(flat_skew_tmp)
      fourth <- sum(a) * (flat_skew_tmp %*% t(flat_skew))

      sigma_j <- n_d / (n * n_star) * (first - second - third + fourth)

      if (j < o) { # force trace to be n_d for all sigmas except last
        sigma_j <- sigma_j / sum(diag(sigma_j)) * n_d
      }

      new_sigmas[[j]] <- sigma_j
    }

    # update all parameters
    mu <- new_mu
    skew <- new_skew
    sigmas <- new_sigmas

    lambda <- new_lambda
    omega <- new_omega

    # Step 5: Check convergence

    total_loglik <- 0

    for(i in 1:n) {
      total_loglik <- total_loglik +
        dtgenhyper(data[[i]], mu, skew, sigmas, lambda, omega, log = TRUE)
    }

    logliks[t] <- total_loglik

    if(t >= 3) {
      ll_rel <- abs(logliks[t] - logliks[t - 1]) / (abs(logliks[t - 1]) + 1e-8)

      if (ll_rel < tol) {
        if (!quiet) message("Converged at iteration ", t)
        break
      }
    }

    if (t %% 50 == 0 & !quiet) {
      cat(sprintf("Iteration %d: criterion = %.3e\n", t, ll_rel))
    }
  }

  if(t == max_iter) message("Reached max iter ", max_iter)

  list(mu = mu, skew = skew, sigmas = sigmas,
       lambda = lambda, omega = omega,
       Ew = a, Einvw = b, Elogw = c)
}
