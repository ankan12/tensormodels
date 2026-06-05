#' tensor_mle_skewt
#'
#' Computes the EM estimation algorithm for the skewt distribution.
#' @return A list containing the estimated mean array, skew array, list of
#'   covariance matrices, and nu.
#'
#' @noRd
tensor_mle_skewt <- function(data, max_iter = 1e3, tol = 1e-6,
                             quiet = TRUE, restrict) {
  # get dim of input
  n <- length(data)
  dims <- dim(data[[1]])
  o <- length(dims)
  n_star <- prod(dims)

  # Step 1: Initialize vals

  # different params based on model
  nu <- 20

  flat_draws <- simplify2array(data)

  mean_draws <- apply(flat_draws, 1:o, mean)
  median_draws <- apply(flat_draws, 1:o, median)
  # res_normal <- tensor_mle(data, model = "normal")
  #
  # mu <- res_normal$mu
  # sigmas <- res_normal$sigmas
  #
  # inv_sigma_start <- lapply(sigmas, invert_safe)
  #
  # precision_resids <- lapply(data, function(x) {
  #   centered <- x - mu
  #
  #   for (d in seq_along(inv_sigma_start)) {
  #     centered <- n_prod(centered, inv_sigma_start[[d]], d)
  #   }
  #
  #   centered
  # })
  #
  # precision_resids <- simplify2array(precision_resids)

  # E[X] = M + nu/(nu-2) * skew
  skew <- (nu-2)/nu * (mean_draws - median_draws)

  logliks <- rep(NA_real_, max_iter)
  pen_logliks <- rep(NA_real_, max_iter)

  sigmas <- lapply(dims, diag)
  mu <- median_draws

  for(k in 1:o) {
    tot_sum <- 0

    for(i in 1:n) {
        curr_unfold <- matricization(data[[i]] - mean_draws, k)

        tot_sum <- tot_sum + (curr_unfold %*% t(curr_unfold))
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

    delta_vals <- delta_vals + nu
    param_vals <- -(nu + n_star) / 2

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

    #replace NaN or Inf values
    # b[!is.finite(b)] <- 0
    # c[!is.finite(c)] <- 0

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

    update_nu <- function(nu, b, c, n) {
      log(nu / 2) + 1 - digamma(nu / 2) - mean(b + c)
    }

    new_nu <- uniroot(
      update_nu,
      interval = c(1e-3, 1e3),
      b = b,
      c = c,
      n = n)$root

    #if(new_nu < 4) new_nu <- 4

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

      new_sigmas[[j]] <- sigma_j
    }

    scale_prod <- 1

    for (j in 1:(o-1)) { # force trace to be n_d for all sigmas except last
      curr_sigma <- new_sigmas[[j]]
      n_d <- dims[j]

      tr_j <- sum(diag(curr_sigma))

      scale_curr <- tr_j / n_d
      scale_prod <- scale_prod * scale_curr

      curr_sigma <- curr_sigma / scale_curr
      new_sigmas[[j]] <- curr_sigma
    }

    new_sigmas[[o]] <- new_sigmas[[o]] * scale_prod

    # update all parameters
    mu <- new_mu
    skew <- new_skew
    sigmas <- new_sigmas
    nu <- new_nu

    # Step 5: Check convergence

    total_loglik <- 0

    for(i in 1:n) {
      total_loglik <- total_loglik +
        dtskewt(data[[i]], mu, skew, sigmas, nu, log = TRUE)
    }

    logliks[t] <- total_loglik

    if(t >= 3) {
      ll_rel <- abs(logliks[t] - logliks[t - 1]) /
                (abs(logliks[t - 1]) + 1e-8)

      if(is.finite(ll_rel) && ll_rel < tol) {
        if(!quiet) message("Converged at iteration ", t)
        break
      }
    }
    if (t %% 50 == 0 & !quiet) {
      cat(sprintf("Iteration %d: criterion = %.3e\n", t, ll_rel))
    }
  }

  if(t == max_iter) message("Reached max iter ", max_iter)

  k <- n_star + sum((dims * (dims+1))/2) - (o - 1) + 1

  list(mu = mu, skew = skew, sigmas = sigmas, nu = nu,
       Ew = a, Einvw = b, Elogw = c,
       loglik = logliks[t],
       k = k,
       AIC = 2 * k - 2 * logliks[t],
       BIC = k * log(n) - 2 * logliks[t])
}
