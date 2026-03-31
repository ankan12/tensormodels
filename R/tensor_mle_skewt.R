#' tensor_mle_skewt
#'
#' Computes the EM estimation algorithm for the skewt distribution.
#' @return A list containing the estimated mean array, skew array, list of
#'   covariance matrices, and nu.
#'
#' @noRd
tensor_mle_skewt <- function(data, max_iter, tol,
                             quiet = TRUE, restrict) {

  # get dim of input
  n <- length(data)
  dims <- dim(data[[1]])
  o <- length(dims)
  n_star <- prod(dims)

  # Step 1: Initialize vals
  mu <- simplify2array(data) |> apply(1:o, mean)

  skew <- array(rnorm(n = n_star), dim = dims)

  logliks <- rep(0, max_iter)

  # s1 <- matrix(rnorm(4), nrow = 2); s2 <- diag(3); s3 <- matrix(rnorm(16), nrow = 4)
  # s1 <- crossprod(s1, s1); s3 <- crossprod(s3, s3)
  #
  est_sigmas <- lapply(dims, diag)

  for(k in 1:o) {
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
  nu <- 4

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

    new_sigmas <- est_sigmas

    # update params based on model
    update_nu <- function(nu, b, c, n) {
      log(nu / 2) + 1 - digamma(nu / 2) - 1 / n * sum(b + c)
    }

    new_nu <- uniroot(
      update_nu,
      interval = c(1e-3, 1e3),
      b = b,
      c = c,
      n = n)$root

    if(new_nu < 2) new_nu <- 2

    scale_prod <- 1

    for (j in 1:o) {
      first <- 0
      second <- 0
      third <- 0
      fourth <- 0

      n_d <- dims[j]

      # covar estimate for mode j
      sigma_j <- matrix(0, nrow = n_d, ncol = n_d)

      # other modes
      other_modes <- (1:o)[-j]

      # Build the whitening operator
      inv_others <- diag(1)

      for (d in rev(other_modes)) {
        inv_others <- kronecker(inv_others, chol2inv(chol(new_sigmas[[d]])))
      }

      A_i <- matricization(new_skew, j)

      for (i in 1:n) {
        Xi_centered <- matricization(data[[i]] - new_mu, j)

        first <- first + b[i] * (Xi_centered %*% inv_others %*% t(Xi_centered))

        second <- second + A_i %*% inv_others %*% t(Xi_centered)

        third <- third + Xi_centered %*% inv_others %*% t(A_i)

        fourth <- fourth + a[i] * A_i %*% inv_others %*% t(A_i)
      }

      sigma_j <- n_d/(n * n_star) * (first - second - third + fourth)

      scale_prod <- 1

      if (j < o) {
        c_d <- sigma_j[1, 1]

        if (!is.finite(c_d) || c_d == 0) c_d <- 1

        scale_prod <- scale_prod * c_d

        sigma_j <- sigma_j / c_d
      }

      new_sigmas[[j]] <- sigma_j
    }

    new_sigmas[[o]] <- new_sigmas[[o]] * scale_prod

    # update all parameters
    mu <- new_mu
    skew <- new_skew
    est_sigmas <- new_sigmas
    nu <- new_nu

    # Step 5: Check convergence

    total_loglik <- 0

    for(i in 1:n) {
      total_loglik <- total_loglik +
        dtskewt(data[[i]], mu, skew, est_sigmas, nu, log = TRUE)
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
    }
    if (t %% 50 == 0 & !quiet) {
      cat(sprintf("Iteration %d: criterion based on Aitken = %.3e\n", t, converge))
    }
  }

  if(t == max_iter) message("Reached max iter ", max_iter)

  k <- n_star + sum((dims * (dims+1))/2) - (o - 1) + 1

  list(mu = mu, skew = skew, sigmas = est_sigmas, nu = nu,
       Ew = a, Einvw = b, Elogw = c,
       loglik = logliks[t], BIC = k * log(n) - 2 * logliks[t])
}
