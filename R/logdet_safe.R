.logdet_safe <- function(S, eps = 1e-8, max_iter = 10, warn = TRUE) {
  p <- nrow(S)
  S <- (S + t(S)) / 2

  sc <- mean(diag(S))
  if (!is.finite(sc) || sc <= 0) {
    sc <- max(abs(S))
  }
  if (!is.finite(sc) || sc <= 0) {
    sc <- 1
  }

  threshold <- eps * sc
  eig_vals <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  min_eig <- min(eig_vals)
  regularized <- FALSE

  if (!is.finite(min_eig) || min_eig < threshold) {
    ridge <- if (is.finite(min_eig)) threshold - min_eig else threshold
    S <- S + diag(ridge, p)
    regularized <- TRUE
  }

  U <- try(chol(S), silent = TRUE)

  if (inherits(U, "try-error")) {
    regularized <- TRUE

    for (k in 0:max_iter) {
      ridge <- eps * (10^k) * sc
      U <- try(chol(S + diag(ridge, p)), silent = TRUE)

      if (!inherits(U, "try-error")) {
        break
      }
    }
  }

  if (inherits(U, "try-error")) {
    stop("Matrix not SPD (even after ridge); check covariance estimation.")
  }

  if (regularized && warn) {
    warning(
      "Covariance matrix is singular or nearly singular; adding ridge regularization before computing the log determinant.",
      call. = FALSE
    )
  }

  2 * sum(log(diag(U)))
}
