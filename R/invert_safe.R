#' invert_safe
#'
#' Inverts a matrix safely. If singular, adds a ridge relative to scale of entries.
#' @return The inverse of the inputted matrix.
#'
#' @noRd
invert_safe <- function(S, eps = 1e-10, max_iter = 6) {
  p <- nrow(S)
  sc <- mean(diag(S))
  if (!is.finite(sc) || sc <= 0) sc <- 1

  for (k in 0:max_iter) {
    lam <- eps * (10^k) * sc
    U <- try(chol(S + diag(lam, p)), silent = TRUE)
    if (!inherits(U, "try-error")) return(chol2inv(U))
  }
  stop("Matrix not SPD (even after ridge); check covariance estimation.")
}
