#' em_est
#'
#' Estimates the mean array and covariance matrices from an array of tensor variate normal draws using an expectation-maximization algorithm.
#'
#' @param draws An array containing the draws, where the first mode represents each draw.
#' @param max_iter A max number of iterations to try to get covariance matrices that converge.
#' @param tol A tolerance level to define the convergence of matrices.
#'
#' @return A list containing the estimated mean array and the list of covariance matrices.
#'
#' @examples
#' s1 <- matrix(rnorm(4), nrow = 2)
#' s1 <- crossprod(s1, s1)
#' s2 <- matrix(rnorm(9), nrow = 3)
#' s2 <- crossprod(s2, s2)
#' s3 <- matrix(rnorm(16), nrow = 4)
#' s3 <- crossprod(s3, s3)
#' rand_mu <- array(rnorm(24), dim = c(2, 3, 4))
#' rtskewt_draws <- rtskewt(n = 1e3, mu = rand_mu, skew = 0.5, sigmas = list(s1, s2, s3))
#' @export
em_est <- function(draws, max_iter = 1000, tol = 1e-6, method = "skewt") {
  # get dim of input
  dims <- dim(draws)[-1]
  num_dim  <- length(dims)
  n <- dim(draws)[1]
  n_star <- prod(dims)

  # Step 1: Initialize vals
  mu <- apply(X = draws, MARGIN = 2:(num_dim+1), FUN = mean)

  skew <- array(1, dim = dims)
  sigmas <- lapply(dims, diag)

  nu <- 20

  fro_rel <- function(Anew, Aold, eps = 1e-12) {
    num <- sqrt(sum((Anew - Aold)^2))
    den <- sqrt(sum(Aold^2)) + eps
    num / den
  }

  max_rel_change <- function(mu_new, mu_old, skew_new, skew_old, sig_new, sig_old, nu_new, nu_old) {
    m_mu   <- fro_rel(mu_new,   mu_old)
    m_skew <- fro_rel(skew_new, skew_old)
    m_sig  <- max(sapply(seq_along(sig_new), function(k) fro_rel(sig_new[[k]], sig_old[[k]])))
    m_nu   <- abs(nu_new - nu_old) / (abs(nu_old) + 1e-12)
    max(c(m_mu, m_skew, m_sig, m_nu))
  }

  for(t in 1:max_iter) {
    # Step 2: Update a, b, c depending on expected values
    skew_compute <- skew

    for(d in seq_along(sigmas)) {
      skew_compute <- n_prod(skew_compute, solve(sigmas[[d]]), d)
    }

    rho <- sum(skew * skew_compute)

    delta_vals <- rep(0, n)

    centered <- draws - array(mu, dim = dim(draws))

    for(i in 1:n) {

      center_draw <- centered[i, , , ]

      centered_compute <- center_draw

      for(d in seq_along(sigmas)) {
        centered_compute <- n_prod(centered_compute, solve(sigmas[[d]]), d)
      }
      delta_vals[i] <- sum(center_draw * centered_compute)
    }

    delta_vals <- delta_vals + nu

    param_vals <- -(nu + n_star)/2

    k_lambda_1 <- besselK(x = sqrt(rho * delta_vals),
                          nu = param_vals + 1,
                          expon.scaled = TRUE)
    k_lambda <- besselK(x = sqrt(rho * delta_vals),
                        nu = param_vals,
                        expon.scaled = TRUE)

     a <- sqrt(delta_vals/rho) * k_lambda_1/k_lambda

     b <- sqrt(rho/delta_vals) * k_lambda_1/k_lambda -
           (2 * param_vals)/delta_vals

     eps = 1e-3

     K_plus <- besselK(x = sqrt(rho * delta_vals),
                       nu = param_vals + eps,
                       expon.scaled = TRUE)
     K_minus <- besselK(x = sqrt(rho * delta_vals),
                        nu = param_vals - eps,
                        expon.scaled = TRUE)
     dK_dlambda <- (K_plus - K_minus) / (2 * eps)

     c <- 0.5 * log(sqrt(delta_vals/rho)) + 1/k_lambda * dK_dlambda

     # ggplot() +
     #   geom_function(fun = update_nu, args = list(b = b, c = c, n = n)) +
     #   geom_hline(yintercept = 0, color = "red", lty = 2) +
     #   xlim(0, 5) +
     #   theme_minimal()

     # Step 3: Update mu, skew, params
     weight_mean <- mean(a) * b - 1
     weight_skew <- mean(b) - b
     num_mean <- 0
     num_skew <- 0

     for (i in 1:n) {
       num_mean <- num_mean + weight_mean[i] * draws[i, , ,]
       num_skew <- num_skew + weight_skew[i] * draws[i, , ,]
     }

     den_mean <- sum(mean(a) * b) - n

     new_mu <- num_mean/den_mean

     den_skew <- sum(a * mean(b)) - n

     new_skew <- num_skew/den_skew

     new_sigmas <- vector("list", num_dim)

     update_nu <- function(nu, b, c, n) {
       log(nu/2) + 1 - digamma(nu/2) - 1/n * sum(b + c)
     }

     new_nu <- uniroot(update_nu, interval = c(1e-3, 1e4), b = b, c = c, n = n)$root

     # Step 4: Update sigmas
     for (j in 1:num_dim) {
       first <- 0; second <- 0; third <- 0; fourth <- 0

       n_d <- dims[j]

       # covar estimate for mode j
       sigma_j <- matrix(0, nrow = n_d, ncol = n_d)

       # other modes
       other_modes <- setdiff(seq_len(num_dim), j)

       # Build the whitening operator (⊗_{d≠j} Δ_d^{-1/2})
       inv_others <- diag(1)

       for(d in rev(other_modes)) {
         inv_others <- kronecker(inv_others, solve(sigmas[[d]]))
       }

       A_i <- matricization(new_skew, j)

       for(i in 1:n) {
          Xi_centered <- matricization(draws[i , , ,] - new_mu, j)

          first <- first + b[i] * (Xi_centered %*% inv_others %*% t(Xi_centered))

          second <- second + A_i %*% inv_others %*% t(Xi_centered)

          third <- third + Xi_centered %*% inv_others %*% t(A_i)

          fourth <- fourth + a[i] * A_i %*% inv_others %*% t(A_i)
       }

       sigma_j <- n_d/(n * n_star) * (first - second - third + fourth)
       sigma_j <- sigma_j / mean(diag(sigma_j))
       new_sigmas[[j]] <- sigma_j
     }

   # Step 5: Check convergence
   max_change <- max_rel_change(new_mu, mu, new_skew, skew, new_sigmas, sigmas, new_nu, nu)

   if(t %% 50 == 0) cat(sprintf("Iteration %d: max relative change = %.3e\n", t, max_change))

   if (max_change < tol) {
     message("Converged at iteration ", t)
     break
   }

   mu <- new_mu
   skew <- new_skew
   sigmas <- new_sigmas
   nu <- new_nu
  }
  list(mu = mu, skew = skew, sigmas = sigmas, nu = nu)
}

#
# ggplot() +
#   geom_function(fun = dinvgamma, args = list(shape = 0.3/2, rate = 0.3/2))

