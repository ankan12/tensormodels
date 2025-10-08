#' nm-mode product
#'
#' Compute the nm-mode product of a tensor with a tensor
#'
#' @param tensorA An array representing the tensor A.
#' @param tensorB An array representing tensor B.
#' @param n An integer specifying which mode of tensor A to multiply across.
#' @param m An integer specifying which mode of tensor B to multiply across.
#'
#' @return An array: the result of the nm-mode product between A and B.
#'
#' @examples
#' a <- array(1:3, dim = c(3, 1, 1))
#' b <- matrix(4:9, nrow = 2, ncol = 3)
#' n_mode_prod(a, b, 1)
#' @export
nm_prod <- function(tensorA, tensorB, n, m) {
  mn_mode_prod_cpp(tensorA, tensorB, n, m)
}
