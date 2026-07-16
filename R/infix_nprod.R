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

#' n-mode product with an explicit mode specification
#'
#' The right-hand side is written as a two-sided formula, with the matrix (or
#' tensor) on the left and the mode on the right. For example,
#' `A %xn% (B ~ 2)` computes the mode-2 product of `A` and `B`.
#'
#' @param tensor A tensor / array.
#' @param specification A two-sided formula of the form `mat ~ mode`.
#'
#' @return An array, or a `tensor` when `tensor` is a tensor object.
#' @export
`%xn%` <- function(tensor, specification) {
  if (!inherits(specification, "formula") || length(specification) != 3L) {
    stop(
      "The right-hand side must be a formula of the form `mat ~ mode`.",
      call. = FALSE
    )
  }

  formula_env <- environment(specification)
  if (is.null(formula_env)) formula_env <- parent.frame()
  mat <- eval(specification[[2L]], envir = formula_env)
  mode <- eval(specification[[3L]], envir = formula_env)

  if (length(mode) != 1L || is.na(mode) || !is.numeric(mode) ||
      !is.finite(mode) || mode != floor(mode) ||
      mode > .Machine$integer.max || mode < 1) {
    stop("`mode` must be a positive integer.", call. = FALSE)
  }

  tensor_is_tensor <- inherits(tensor, "tensor")
  tensor_input <- if (tensor_is_tensor) unclass(tensor) else tensor
  mat_input <- if (inherits(mat, "tensor")) unclass(mat) else mat
  result <- n_prod(tensor_input, mat_input, n = as.integer(mode))

  if (tensor_is_tensor) tensor(result) else result
}
