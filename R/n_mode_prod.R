#' n-mode product
#'
#' Compute the n-mode product of a tensor with a matrix.
#'
#' @param tensor An array representing the tensor.
#' @param mat A matrix of size \eqn{m \times n} to multiply with.
#' @param n  An integer specifying the mode to multiply across.
#'
#' @return An array: the result of the n-mode product.
#'
#' @examples
#' a <- array(1:3, dim = c(3, 1, 1))
#' b <- matrix(4:9, nrow = 2, ncol = 3)
#' n_mode_prod(a, b, 1)
#' @export
n_mode_prod <- function(tensor, mat, n) {
  n_mode_prod_cpp(tensor, mat, n)
}
