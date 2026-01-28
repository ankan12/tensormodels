#' frob_norm
#'
#' Computes the Frobenius norm of a tensor.
#'
#' @param A An array containing the tensor.
#'
#' @return The Frobenius norm.
#' @details The norm of a tensor is the square root of the sum of the
#' squares of all its elements.
#' \deqn{||\mathcal{X}|| = \sqrt{\sum_{i_{1}=1}^{I_{1}} \sum_{i_{2}=1}^{I_{2}} \dots \sum_{i_{N}=1}^{I_{N}} x_{i_{1}, i_{2}, \dots, i_{N}}^{2}}}
#'
#' @examples
#' A <- array(1:24, dim = c(2, 3, 4))
#' frob_norm(A)
#'
#' @export
frob_norm <- function(A) {
  sqrt(sum(as.vector(A)^2))
}
