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
#' @export

em_est <- function(draws, max_iter = 1000, tol = 1e-6, quiet = TRUE, model) {

  if(!model %in% c("skewt", "vargamma", "invgamma", "genhyper")) stop("Not a valid model. Must be skewt, vargamma, invgamma, or genhyper")

  # get dim of input
  dims <- dim(draws)[-1]
  num_dim <- length(dims)
  n <- dim(draws)[1]
  n_star <- prod(dims)

  # Step 1: Initialize vals
  mu <- apply(X = draws, MARGIN = 2:(num_dim+1), FUN = mean)

  skew <- array(rnorm(prod(dims)), dim = dims)
  sigmas <- lapply(dims, diag)

  # different params based on model
  switch(model,
         "skewt" = nu <- 10,
         "vargamma" = gamma <- 10,
         "invgamma" = kappa <- 10,
         "genhyper" = {lambda <- 10; omega <- 10 })

  for(t in 1:max_iter) {
    # Step 2: Update a, b, c depending on expected values
    skew_compute <- skew

    for(d in seq_along(sigmas)) {
      skew_compute <- n_prod(skew_compute, chol2inv(chol(sigmas[[d]])), d)
    }

    rho <- sum(skew * skew_compute)

    delta_vals <- rep(0, n)

    centered <- draws - array(mu, dim = dim(draws))

    for(i in 1:n) {

      center_draw <- centered[i, , , ]

      centered_compute <- center_draw

      for(d in seq_along(sigmas)) {
        centered_compute <- n_prod(centered_compute, chol2inv(chol(sigmas[[d]])), d)
      }
      delta_vals[i] <- sum(center_draw * centered_compute)
    }

    switch(model,
           "skewt" = {
             delta_vals <- delta_vals + nu
             param_vals <- -(nu + n_star)/2},
           "genhyper" = {
             rho <- rho + omega
             delta_vals <- delta_vals + omega
             param_vals <- lambda - n_star/2
            },
           "vargamma" = {
             rho <- rho + 2 * gamma
             param_vals <- gamma - n_star/2
            },
           "invgamma" = {
             rho <- rho + kappa^2
             delta_vals <- delta_vals + 1
             param_vals <- -(1 + n_star)/2
           })

    k_lambda_1 <- besselK(x = sqrt(rho * delta_vals),
                          nu = param_vals + 1, expon.scaled = TRUE)
    k_lambda <- besselK(x = sqrt(rho * delta_vals),
                        nu = param_vals, expon.scaled = TRUE)

    a <- sqrt(delta_vals/rho) * (k_lambda_1/k_lambda)

    b <- sqrt(rho/delta_vals) * (k_lambda_1/k_lambda) -
         (2 * param_vals)/delta_vals

    eps <- 1e-5

    K_plus <- besselK(x = sqrt(rho * delta_vals), nu = param_vals + eps,
                      expon.scaled = TRUE)

    K_minus <- besselK(x = sqrt(rho * delta_vals), nu = param_vals - eps,
                       expon.scaled = TRUE)

    c <- 1/2 * log(delta_vals/rho) + (log(K_plus) - log(K_minus))/(2 * eps)

    # replace NaN or Inf values
    b[!is.finite(b)] <- 0
    c[!is.finite(c)] <- 0

    #browser()
    # Step 3: Update mu, skew, params
    weight_mean <- mean(a) * b - 1
    weight_skew <- mean(b) - b

    num_mean <- 0
    num_skew <- 0

    for (i in 1:n) {
       num_mean <- num_mean + weight_mean[i] * draws[i, , ,]
       num_skew <- num_skew + weight_skew[i] * draws[i, , ,]
    }

    den_mean <- sum(mean(a) * b - n)

    new_mu <- num_mean/den_mean

    den_skew <- sum(a * mean(b) - n)
    new_skew <- num_skew/den_skew

    new_sigmas <- sigmas

    # update params based on model

    update_nu <- function(nu, b, c, n) log(nu/2) + 1 - digamma(nu/2) - 1/n * sum(b+c)
    update_gamma <- function(gamma, a, c) log(gamma) + 1 - digamma(gamma) + mean(c) - mean(a)


    switch(model,
           "skewt" = {
             new_nu <- uniroot(update_nu, interval = c(1e-3, 1e3),
                               b = b, c = c, n = n)$root
            },
           "genhyper" = {
             K_plus <- besselK(x = omega, nu = lambda + eps,
                               expon.scaled = TRUE)

             K_minus <- besselK(x = omega, nu = lambda - eps,
                                expon.scaled = TRUE)

             new_lambda <- mean(c) * lambda * 1/((log(K_plus) - log(K_minus))/(2 * eps))

             R_lambda <-
               besselK(x = omega, nu = new_lambda + 1, expon.scaled = TRUE)/
               besselK(x = omega, nu = new_lambda, expon.scaled = TRUE)

             R_neg_lambda <-
               besselK(x = omega, nu = -new_lambda + 1, expon.scaled = TRUE)/
               besselK(x = omega, nu = -new_lambda, expon.scaled = TRUE)

             first_deriv <- 1/2(R_lambda + R_neg_lambda - (mean(a) + mean(b)))

             second_deriv <- R_lambda^2 - (1 + 2 * new_lambda)/omega * R_lambda - 1 +
                                R_neg_lambda^2 - (1 - 2 * lambda)/omega * R_neg_lambda - 1

             new_omega <- omega - first_deriv/second_deriv
           },
           "vargamma" = {
             new_gamma <- uniroot(update_gamma, interval = c(1e-3, 1e3),
                                  a = a, c = c)$root
           },
           "invgamma" = {
             new_kappa = n/(sum(a))
           })

    scale_prod <- 1

  for (j in 1:num_dim) {
    first <- 0; second <- 0; third <- 0; fourth <- 0

    n_d <- dims[j]

    # covar estimate for mode j
    sigma_j <- matrix(0, nrow = n_d, ncol = n_d)

    # other modes
    other_modes <- (1:num_dim)[-j]

    # Build the whitening operator
    inv_others <- diag(1)

    for(d in rev(other_modes)) {
      inv_others <- kronecker(inv_others, chol2inv(chol(sigmas[[d]])))
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

    scale_prod <- 1

    if(j < num_dim) sigma_j <- sigma_j/(sum(diag(sigma_j))) * n_d

    new_sigmas[[j]] <- sigma_j
  }

   #new_mu <- rand_mu; new_skew <- new_skew; new_sigmas <- list(s1, s2, s3); new_nu <- 20

   # Step 5: Check convergence
   mean_change <- mean_rel_change(new_mu, mu, new_skew, skew, new_sigmas, sigmas, new_nu, nu)

   if(t %% 50 == 0 & !quiet) {
     cat(sprintf("Iteration %d: mean relative change = %.3e\n", t, mean_change))
     #print(list(mu = mu, skew = skew, sigmas = sigmas, nu = nu))
   }

   if (mean_change < tol) {
     message("Converged at iteration ", t)
     break
   }
   # update all parameters
   mu <- new_mu; skew <- new_skew; sigmas <- new_sigmas;

   switch(model,
          "skewt" = nu <- new_nu,
          "genhyper" = {lambda <- new_lambda; omega <- new_omega},
          "vargamma" = gamma <- new_gamma,
          "invgamma" = kappa <- new_kappa)
  }
  switch(model,
         "skewt" = list(mu = mu, skew = skew, sigmas = sigmas, nu = nu),
         "genhyper" = list(mu = mu, skew = skew, sigmas = sigmas, lambda = lambda, omega = omega),
         "vargamma" = list(mu = mu, skew = skew, sigmas = sigmas, gamma = gamma),
         "invgamma" = list(mu = mu, skew = skew, sigmas = sigmas, kappa = kappa))
}

# Less computationally intensive way to compute sigmas

# Step 4: Update sigmas
# n_2 <- dims[[2]]
#
# n3_star <- prod(dims[3:num_dim])
#
# inv_sigma2 <- chol2inv(chol(new_sigmas[[2]]))
#
# sigma_1 <- 0
# sigma_2 <- 0
# curr_sigma <- 0
#
# for(i in 1:n) {
#   Xi_centered <- matricization(draws[i , , ,] - new_mu, 1)
#
#   for(j in 1:n3_star) {
#     # e_j <- rep(0, nrow = dims[[3]])
#     # e_j[j] <- 1
#
#     #e_j <- rep(0, n3_star)
#     #e_j[j] <- 1
#     #e_j <- matrix(e_j, ncol = 1)
#     #e_j <- diag(n_2)[, j, drop = FALSE]
#
#     A1j <- diag(n_2)
#     X1ij <- diag(n_2)
#    #A1j <- kronecker(diag(n_2), t(e_j)) %*% matricization(new_skew, 1)
#     #X1ij <- kronecker(diag(n_2), t(e_j))
#
#    # multiply by other covar past 3
#    for(d in 3:num_dim) {
#       A1j <- kronecker(A1j, chol2inv(chol(sigmas[[d]])))
#       X1ij <- kronecker(X1ij, chol2inv(chol(sigmas[[d]])))
#    }
#
#    e_j <- matrix(0, nrow = n3_star, ncol = 1)#nrow = dims[[3]], ncol = 1)
#    e_j[j,1] <- 1
#
#    A1j <- (kronecker(diag(n_2), t(e_j)) %*% A1j) %*% t(matricization(new_skew, 1))
#    X1ij <- (kronecker(diag(n_2), t(e_j)) %*% X1ij) %*% t(matricization(Xi_centered, 1))
#
#    #A1j <- kronecker(A1j, e_j) %*% matricization(new_skew, 1)
#    #X1ij <- kronecker(diag(n_2), t(e_j)) %*% Xi_centered
#
#    #A1j <- A1j %*%
#    #X1ij <- X1ij %*% Xi_centered
#
#    sigma_1 <- sigma_1 + -t(A1j) %*% inv_sigma2 %*% X1ij -
#                          t(X1ij) %*% inv_sigma2 %*% A1j +
#                          b[i] * t(X1ij) %*% inv_sigma2 %*% X1ij +
#                          a[i] * t(A1j) %*% inv_sigma2 %*% A1j
#    #A_i <- matricization(new_skew, j)
#   }
# }
#   sigma_1 <- sigma_1 * dims[[1]]/(n * n_star)
#
#   new_sigmas[[1]] <- sigma_1/sum(diag(sigma_1))
#
#   # Solve for sigma2
#
#   for(i in 1:n) {
#     Xi_centered <- matricization(draws[i , , ,] - new_mu, 1)
#
#     for(j in 1:n3_star) {
#       # e_j <- rep(0, n3_star)
#       # e_j[j] <- 1
#       # e_j <- matrix(e_j, ncol = 1)
#       # #e_j <- diag(n_2)[, j, drop = FALSE]
#       #
#       # A1j <- kronecker(diag(n_2), t(e_j))
#       # X1ij <- kronecker(diag(n_2), t(e_j))
#
#       A1j <- diag(n_2)
#       X1ij <- diag(n_2)
#
#       # multiply by other covar past 3
#       for(d in 3:num_dim) {
#         A1j <- kronecker(A1j, chol2inv(chol(sigmas[[d]])))
#         X1ij <- kronecker(X1ij, chol2inv(chol(sigmas[[d]])))
#       }
#       # A1j <- A1j %*% matricization(new_skew, 1)
#       # X1ij <- X1ij %*% Xi_centered
#
#       e_j <- matrix(0, nrow = n3_star, ncol = 1)
#       e_j[j,1] <- 1
#
#       A1j <- (kronecker(diag(n_2), t(e_j)) %*% A1j) %*% t(matricization(new_skew, 1))
#       X1ij <- (kronecker(diag(n_2), t(e_j)) %*% X1ij) %*% t(matricization(Xi_centered, 1))
#
#       inv_sigma1 <- chol2inv(chol(sigmas[[1]]))
#
#       sigma_2 <- sigma_2 + -(X1ij) %*% inv_sigma1 %*% t(A1j) -
#                             A1j %*% inv_sigma1 %*% t(X1ij) +
#                             b[i] * X1ij %*% inv_sigma1 %*% t(X1ij) +
#                             a[i] * A1j %*% inv_sigma1 %*% t(A1j)
#
#       # sigma_1 <- sigma_1 + -t(A1j) %*% chol2inv(chol(new_sigmas[[2]])) * X1ij -
#       #   t(X1ij) %*% chol2inv(chol(new_sigmas[[2]])) * A1j +
#       #   b[i] * t(X1ij) %*% chol2inv(chol(new_sigmas[[2]])) %*% X1ij +
#       #   a[i] * t(A1j) %*% chol2inv(chol(new_sigmas[[2]])) %*% A1j
#       #A_i <- matricization(new_skew, j)
#     }
#   }
#   sigma_2 <- sigma_2 * dims[[2]]/(n * n_star)
#
#   new_sigmas[[2]] <- sigma_2/sum(diag(sigma_2))
#
#   # Solve for all other sigmas
#
#   for(l in 3:num_dims) {
#     curr_sigma <- 0
#
#     n2_star_Dl <- prod(dims[-c(2, l)])
#     n_l <- dims[[l]]
#     for(i in 1:N) {
#       for(j in 1:n2_star_Dl) {
#         A1jl2 <- diag(n_l)
#         X1ijl2 <- diag(n_l)
#
#         # multiply by other covar past 2
#         for(d in 2:num_dim) {
#           if(d != l) { # skip the lth covar
#             A1jl2 <- kronecker(A1jl2, chol2inv(chol(sigmas[[d]]))
#             X1ijl2 <- kronecker(X1ijl2, chol2inv(chol(sigmas[[d]]))
#           }
#         }
#         # A1j <- A1j %*% matricization(new_skew, 1)
#         # X1ij <- X1ij %*% Xi_centered
#
#         e_j <- matrix(0, nrow = 3, ncol = 1)
#         e_j[j,1] <- 1
#
#         A1jl2 <- (kronecker(diag(n_l), t(e_j)) %*% A1jl2) %*% t(matricization(new_skew, 1))
#         X1ijl2 <- (kronecker(diag(n_l), t(e_j)) %*% X1ijl2) %*% t(matricization(Xi_centered, 1))
#
#         inv_sigma1 <- chol2inv(chol(sigmas[[1]]))
#
#         curr_sigma <- curr_sigma +
#                       -(X1ijl2) %*% inv_sigma1 %*% t(A1jl2) -
#                        A1jl2 %*% inv_sigma1 %*% t(X1ijl2) +
#                        b[i] * X1ijl2 %*% inv_sigma1 %*% t(X1ijl2) +
#                        a[i] * A1jl2 %*% inv_sigma1 %*% t(A1jl2)
#       }
#     }
#
#     curr_sigma <- curr_sigma * dims[[l]]/(n * n_star)
#
#     new_sigmas[[l]] <- curr_sigma/sum(diag(curr_sigma))
#   }
