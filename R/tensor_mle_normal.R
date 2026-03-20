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
tensor_mle_normal <- function(data, max_iter = 1000, tol = 1e-6,
                              quiet = TRUE, restrict = NULL) {
  #get dim of input
  n <- length(data)
  dims <- dim(data[[1]])
  num_dim <- length(dims)

  if (length(restrict) >= num_dim) {
    stop("Invalid restriction: You must restrict at most number of dims - 1 scale parameters.")
  }

  if(num_dim == 1 && dims == 1) {
    vector_data <- simplify2array(data)
    mu <- mean(vector_data)
    sigma <- sd(vector_data)
    return(list(mu = mu, sigmas = list(sigma)))
  }

  if(num_dim == 1) {
    array_data <- simplify2array(data) |> aperm(c(2, 1))
    sigma <- cov(array_data) * (n - 1)/n
    return(list(mu = mu, sigmas = list(sigma)))
  }

  mu <- simplify2array(data) |> apply(1:num_dim, mean)

  #intialize sigmas as identity matrices
  est_sigmas <- lapply(dims, diag)

  logliks <- rep(0, max_iter)

  for (t in 1:max_iter) {
    # store previous iteration
    old_sigmas <- lapply(est_sigmas, identity)

    for (k in 1:num_dim) {
      #inverses of all sigma
      inv_sigma <- lapply(est_sigmas, invert_safe)

      # dim for kth mode
      d_k    <- dims[k]
      d_negk <- prod(dims[-k])
      s_k    <- matrix(0, d_k, d_k)

      # accumlate mode-k covar
      for (draw in 1:n) {
        Xi <- data[[draw]] - mu
        Xik <- matricization(Xi, k)

        for(index_kroneck in (1:num_dim)[-k]) {
          Xi <- n_prod(Xi, inv_sigma[[index_kroneck]], index_kroneck)
        }

        s_k <- s_k + matricization(Xi, k) %*% t(Xik)
      }

      if (k %in% restrict) { # restrict to be product of identity
        s_k <- s_k / (n * d_negk)
        curr_sigma <- sum(diag(s_k)) / d_k
        est_sigmas[[k]] <- diag(d_k) * curr_sigma
        next
      }
      # update sigma_k
      est_sigmas[[k]] <- s_k / (n * d_negk)
    }

    scale_prod <- 1

    for (d in 1:(num_dim - 1)) {
      if (d %in% restrict) next  # don't rescale constrained factors

      c_d <- est_sigmas[[d]][1, 1]

      if (!is.finite(c_d) || c_d == 0) c_d <- 1

      est_sigmas[[d]] <- est_sigmas[[d]] / c_d
      scale_prod <- scale_prod * c_d
    }

    est_sigmas[[num_dim]] <- est_sigmas[[num_dim]] * scale_prod

    # Step 5: Check convergence

    total_loglik <- 0

    for(i in 1:n) {
      total_loglik <- total_loglik +
                      dtnorm(data[[i]], mu, est_sigmas, log = TRUE)
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

  list(mu = mu, sigmas = est_sigmas)
}
