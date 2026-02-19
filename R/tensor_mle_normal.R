#' tensor_mle_normal
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
#' fourth_tensor <- rtnorm(n = 1e3, mu = fourth_mu, sigmas = fourth_sigma)
#' est_fourth <- mle_est(fourth_tensor)
#' (fourth_mu - est_fourth$mu)^2 |> mean()
#' fourth_scaled <- lapply(fourth_sigma, function(S) S * (nrow(S) / sum(diag(S))))
#' (mse_each <- mapply(function(est, true) mean((est - true)^2), est_fourth$sigmas, fourth_scaled))
#' @export
tensor_mle_normal <- function(data, max_iter = 1000, tol = 1e-6, quiet = TRUE) {
  #get dim of input
  all_dims <- dim(data)[-1]
  num_dim  <- length(all_dims)
  n <- dim(data)[1]

  if(num_dim == 1) {
    mu <- apply(data, 2, mean)
    sigma <- cov(data) * (n - 1)/n
    #sigma <- sigma * (nrow(sigma) / sum(diag(sigma)))
    return(list(mu = mu, sigmas = list(sigma)))
  }

  #center input data
  mean_array <- apply(data, MARGIN = 2:(num_dim+1), FUN = mean)
  centered_X <- sweep(data, MARGIN = 2:(num_dim+1), STATS = mean_array, FUN = "-")

  #intialize as identity matrices
  est_sigmas <- lapply(all_dims, diag)

  logliks <- rep(0, max_iter)

  for (t in 1:max_iter) {
    # store previous iteration
    old_sigmas <- lapply(est_sigmas, identity)

    for (k in 1:num_dim) {
      #inverses of all sigma
      inv_sigma <- lapply(est_sigmas, solve)

      #exclude kth mode
      inv_sigma_except_k <- inv_sigma[-k]

      # Kronecker product of inverses except k
      kron_except_k <- Reduce(function(A,B) kronecker(B, A), inv_sigma_except_k)

      # dim for kth mode
      d_k    <- all_dims[k]
      d_negk <- prod(all_dims[-k])
      s_k    <- matrix(0, d_k, d_k)

      # accumlate mode-k covar
      for (draw in 1:n) {
        idx_list <- c(list(draw), rep(list(bquote()), num_dim))
        Xi <- do.call("[", c(list(centered_X), idx_list, list(drop = TRUE)))
        Xik <- tensormodels::matricization(Xi, k)
        s_k <- s_k + Xik %*% kron_except_k %*% t(Xik)
      }

      # update sigma_k
      est_sigmas[[k]] <- s_k / (n * d_negk)

      # normalize sigmas
      est_sigmas[[k]] <- est_sigmas[[k]] / sum(diag(est_sigmas[[k]]))
    }

    # Step 5: Check convergence

    logliks[t] <- loglik_normal_observed(data, mean_array, est_sigmas)

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
    # convergence check
    # rel_changes <- mapply(function(new, old)
    #   norm(new - old, "F") / (norm(old, "F") + 1e-12),
    #   est_sigmas, old_sigmas)
    #
    # max_change <- max(rel_changes)
    #
    # if (t %% 50 == 0 & !quiet) {
    #   cat(sprintf(
    #     "Iteration %d: mean relative change = %.3e\n",
    #     t,
    #     mean_change
    #   ))
    # }
    #
    # if (max_change < tol) {
    #   if(!quiet) message("Converged at iteration ", t)
    #   break
    # }
  }

  # normalize covariance matrices with trace
  list(mu = mean_array,
       sigmas = lapply(est_sigmas, function(S) S * (nrow(S) / sum(diag(S)))))
}
