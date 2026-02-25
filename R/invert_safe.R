#' invert_safe
#'
#' Inverts a matrix safely. If singular, adds a ridge relative to scale of entries.
#' @return The inverse of the inputted matrix.
#'
#' @noRd
invert_safe <- function(S, eps = 1e-10, max_iter = 10) {
  chol2inv(chol_safe(S, eps = eps, max_iter = max_iter))
}
