#' mle_est
#'
#' Estimates the mean array and covariance matrices from an array of tensor variate normal draws.
#'
#' @param data An array containing the draws, where the first mode represents each draw.
#' @param max_iter A max number of iterations to try to get covariance matrices that converge.
#' @param tol A tolerance level to define the convergence of matrices.
#'
#' @return A list containing the estimated mean array and the list of covariance matrices.
#'
#' @examples
#' all_dims <- c(2, 10, 10, 2)
#' fourth_mu = array(rnorm(400), dim = all_dims)
#' fourth_sigma = lapply(seq_along(all_dims), function(k) diag(all_dims[k]))
#' fourth_tensor <- rtnorm(n = 1e3, mu = fourth_mu, listSigmas = fourth_sigma)
#' est_fourth <- mle_est(fourth_tensor)
#' (fourth_mu - est_fourth$mu)^2 |> mean()
#' fourth_scaled <- lapply(fourth_sigma, function(S) S / det(S)^(1/nrow(S)))
#' (mse_each <- mapply(function(est, true) mean((est - true)^2), est_fourth$sigmas, fourth_scaled))
#' @export
mle_est <- function(data, max_iter = 1000, tol = 1e-6) {
  #get dim of input
  all_dims <- dim(data)[-1]
  num_dim  <- length(all_dims)
  n <- dim(data)[1]

  if (num_dim == 1) { #univar case
    n <- length(data)
    mu <- mean(data)
    sigma <- var(data)
    return(list(mu = mu, sigmas = list(matrix(sigma, 1, 1))))
  }

  #center input data
  mean_array <- apply(data, MARGIN = 2:(num_dim+1), FUN = mean)
  centeredX <- sweep(data, MARGIN = 2:(num_dim+1), STATS = mean_array, FUN = "-")

  #intialize as identity matrices
  est_sigmas <- lapply(all_dims, diag)

  for (t in 1:max_iter) {
    old_sigmas <- lapply(est_sigmas, identity)  # store previous iteration

    for (k in 1:num_dim) {
      invSigma <- lapply(est_sigmas, solve) #inverses of all sigma

      invSigma_except_k <- invSigma[-k] #exclude kth mode

      #Kronecker product of inverses except k
      kron_except_k <- Reduce(function(A,B) kronecker(B, A), invSigma_except_k)

      #dim for kth mode
      d_k    <- all_dims[k]
      d_negk <- prod(all_dims[-k])
      s_k    <- matrix(0, d_k, d_k)

      #accumlate mode-k covar
      for (draw in 1:n) {
        idx_list <- c(list(draw), rep(list(bquote()), num_dim))
        Xi <- do.call("[", c(list(centeredX), idx_list, list(drop = TRUE)))
        Xik <- tensormodels::matricization(Xi, k)
        s_k <- s_k + Xik %*% kron_except_k %*% t(Xik)
      }

      #update sigma_k
      est_sigmas[[k]] <- s_k / (n * d_negk)

      #normalize sigmas
      est_sigmas[[k]] <- est_sigmas[[k]] / mean(diag(est_sigmas[[k]]))
    }

    #convergence check
    rel_changes <- mapply(function(new, old)
      norm(new - old, "F") / (norm(old, "F") + 1e-12),
      est_sigmas, old_sigmas)

    max_change <- max(rel_changes)
    cat(sprintf("Iteration %d: max relative change = %.3e\n", t, max_change))

    if (max_change < tol) {
      message("Converged at iteration ", t)
      break
    }
  }
  list(mu = mean_array, sigmas =
         lapply(est_sigmas, function(S) S / det(S)^(1/nrow(S))))
}
