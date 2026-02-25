#' chol_safe
#'
#' Gets the cholesky a matrix safely.
#' If singular, adds a ridge relative to scale of entries.
#' @return The cholesky of the inputted matrix.
#'
#' @noRd
chol_safe <- function(S, eps = 1e-10, max_iter = 10) {
  U <- try(chol(S), silent = TRUE) # try chol on U

  if (!inherits(U, "try-error")) {
    return(U)
  }

  # error, add ridge to prevent singular matrix
  p <- nrow(S)
  sc <- mean(diag(S))
  if (!is.finite(sc) || sc <= 0) sc <- 1

  for (k in 0:max_iter) {
    lam <- eps * (10^k) * sc
    U <- try(chol(S + diag(lam, p)), silent = TRUE)
    if (!inherits(U, "try-error")) return(U)
  }

  stop("Matrix not SPD (even after ridge); check covariance estimation.")
}
