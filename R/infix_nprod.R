#' n-mode product (mode 1)
#'
#' First mode of tensor times columns of matrix.
#'
#' @param tensor A tensor / array.
#' @param mat A matrix with compatible dimensions.
#' @export
`%x1%` <- function(tensor, mat) {
  n_prod(tensor, mat, n = 1)
}

#' n-mode product (mode 2)
#'
#' Second mode of tensor times columns of matrix.
#'
#' @param tensor A tensor / array.
#' @param mat A matrix with compatible dimensions.
#' @export
`%x2%` <- function(tensor, mat) {
  n_prod(tensor, mat, n = 2)
}

#' n-mode product (mode 3)
#'
#' Third mode of tensor times columns of matrix.
#'
#' @param tensor A tensor / array.
#' @param mat A matrix with compatible dimensions.
#' @export
`%x3%` <- function(tensor, mat) {
  n_prod(tensor, mat, n = 3)
}
