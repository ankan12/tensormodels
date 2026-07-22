#' n-mode product
#'
#' @param tensor A tensor object or an array.
#' @param mat A compatible matrix.
#' @param n The mode within each tensor draw.
#' @return A tensor when `tensor` is a tensor object, otherwise an array.
#' @export
n_prod <- function(tensor, mat, n) {
  tensor_is_tensor <- inherits(tensor, "tensor")
  if (inherits(mat, "tensor")) mat <- .tensor_single_draw_array(mat, "mat")
  if (!is.matrix(mat)) stop("`mat` must be a matrix.", call. = FALSE)
  if (length(n) != 1L || is.na(n) || !is.numeric(n) || !is.finite(n) ||
      n != floor(n) || n < 1L) stop("`n` must be a positive integer.", call. = FALSE)
  n <- as.integer(n)
  if (tensor_is_tensor) {
    if (n > length(draw_shape(tensor))) stop("Invalid mode index.", call. = FALSE)
    out <- .Call(`_tensormodels_n_prod`, unclass(tensor), mat, n + 1L)
    return(.new_tensor_array(out))
  }
  .Call(`_tensormodels_n_prod`, tensor, mat, n)
}
