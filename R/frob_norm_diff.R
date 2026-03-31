#' frob_norm_diff
#'
#' Computes the relative difference of Frobenius norms between two matrices.
#'
#' @param A A matrix.
#' @param B A matrix.
#'
#' @return A value representing the relative Frobenius difference between A and B.
#'
#' @examples
#' matA <- matrix(1:4, nrow = 2)
#' matB <- matrix(1:4, nrow = 2) + rnorm(4)
#' frob_norm_diff(matA, matB)
#'
#' @export
frob_norm_diff <- function(A, B) {
  num <- sqrt(sum((A - B)^2)) # Frobenius norm of the difference
  den <- sqrt(sum(B^2)) # Frobenius norm of reference
  num / den # relative Frobenius difference
}
