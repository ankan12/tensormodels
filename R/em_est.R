#' em_est
#'
#' Estimates the mean array and covariance matrices from an array of tensor variate normal draws using an expectation-maximization algorithm.
#'
#' @param data An array containing the draws, where the first mode represents each draw.
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
#' rand_mu <- array(rtnorm(24), dim = c(2, 3, 4))
#' rtskewt_draws <- rtskewt(n = 1e3, mu = rand_mu, sigmas = list(s1, s2, s3))
#' @export
em_est <- function(data, max_iter = 1000, tol = 1e-6, method = "skewt") {
  # get dim of input
  dims <- dim(data)[-1]
  num_dim  <- length(all_dims)
  n <- dim(data)[1]
  n_star <- prod(dims)

  # initialize vals
  mu <- apply(X = data, MARGIN = 2:num_dim, FUN = mean)
  skew <- array(1, dim = dims)
  sigmas <- lapply(dims, diag)
  nu = 1

  a = 1
  b = 1

  #rho is vec(A)^T otimes Sigma^{-1} vec{A}

  while(diff < tol) {
#
#     rho_vals <- skew
#     for (d in seq_along(sigmas)) {
#       rho_vals <- n_prod(rho_vals, solve(chol(sigmas[[d]])), d)
#       centered <- n_prod(centered, solve(chol(sigmas[[d]])), d)
#     }
#     rho <- sum(rho_vals^2)
#

    rho_vals <- skew
    for (d in seq_along(sigmas)) {
      rho_vals <- n_prod(rho_vals, solve(chol(sigmas[[d]])), d)
    }
    rho <- sum(rho_vals^2) # constant rho for all draws

    delta_vals <- rep(1, n)
    centered <- data - array(mu, dim = dim(data))

    for (i in 1:n) {
      center_draw <- centered[i, , , ]
      for (d in seq_along(sigmas)) {
        center_draw <- n_prod(center_draw, solve(chol(sigmas[[d]])), d)
      }
      delta_vals[i] <- sum(center_draw^2) # different delta for each draw
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

   eps = 1e-2

   K_plus <- besselK(x = sqrt(rho * delta_vals), nu = param_vals + 1 + eps)
   K_minus <- besselK(x = sqrt(rho * delta_vals), nu = param_vals + 1 - eps)
   dK_dlambda <- (K_plus - K_minus) / (2 * eps)

   c <- log(sqrt(delta_vals/rho)) + 1/k_lambda * dK_dlambda

   update_nu <- function(nu, b, c, n) {
     log(nu/2) + 1 - digamma(nu/2) - 1/n * sum(b + c)
   }

   nu <- uniroot(update_nu, interval = c(1e-3, 200), b = b, c = c, n = n)$root
  }

#
#   if (num_dim == 1) { #univar case
#     n <- length(data)
#     mu <- mean(data)
#     sigma <- var(data)
#     return(list(mu = mu, sigmas = list(matrix(sigma, 1, 1))))
#   }
}


#
#   #center input data
#   mean_array <- apply(data, MARGIN = 2:(num_dim+1), FUN = mean)
#   centered_X <- sweep(data, MARGIN = 2:(num_dim+1), STATS = mean_array, FUN = "-")
#
#   #intialize as identity matrices
#   est_sigmas <- lapply(all_dims, diag)
#
#   for (t in 1:max_iter) {
#     # store previous iteration
#     old_sigmas <- lapply(est_sigmas, identity)
#
#     for (k in 1:num_dim) {
#       #inverses of all sigma
#       inv_sigma <- lapply(est_sigmas, solve)
#
#       #exclude kth mode
#       inv_sigma_except_k <- inv_sigma[-k]
#
#       # Kronecker product of inverses except k
#       kron_except_k <- Reduce(function(A,B) kronecker(B, A), inv_sigma_except_k)
#
#       # dim for kth mode
#       d_k    <- all_dims[k]
#       d_negk <- prod(all_dims[-k])
#       s_k    <- matrix(0, d_k, d_k)
#
#       # accumlate mode-k covar
#       for (draw in 1:n) {
#         idx_list <- c(list(draw), rep(list(bquote()), num_dim))
#         Xi <- do.call("[", c(list(centered_X), idx_list, list(drop = TRUE)))
#         Xik <- tensormodels::matricization(Xi, k)
#         s_k <- s_k + Xik %*% kron_except_k %*% t(Xik)
#       }
#
#       # update sigma_k
#       est_sigmas[[k]] <- s_k / (n * d_negk)
#
#       # normalize sigmas
#       est_sigmas[[k]] <- est_sigmas[[k]] / mean(diag(est_sigmas[[k]]))
#     }
#
#     # convergence check
#     rel_changes <- mapply(function(new, old)
#       norm(new - old, "F") / (norm(old, "F") + 1e-12),
#       est_sigmas, old_sigmas)
#
#     max_change <- max(rel_changes)
#     cat(sprintf("Iteration %d: max relative change = %.3e\n", t, max_change))
#
#     if (max_change < tol) {
#       message("Converged at iteration ", t)
#       break
#     }
#   }
#
#   # normalize covariance matrices with trace
#   list(mu = mean_array,
#        sigmas = lapply(est_sigmas, function(S) S * (nrow(S) / sum(diag(S)))))
# }
