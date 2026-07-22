#' Subset a tensor
#'
#' Tensor subsetting follows ordinary array rules. For a one-draw tensor, the
#' observation axis is implicit, so `x[i, j, k]` indexes the tensor modes. For
#' multiple draws, the first index is the draw axis, as in the underlying
#' array. Use [pull_draw()] when an explicit one-draw tensor is desired.
#'
#' @param x A tensor object.
#' @param ... Array indices.
#' @param drop Whether singleton dimensions should be dropped.
#' @return The subsetted value, retaining class when the result is an array
#'   with an explicit draw axis.
#' @export
#' @method [ tensor
`[.tensor` <- function(x, ..., drop = TRUE) {
  args <- as.list(substitute(list(...)))[-1L]
  if (length(args) == 0L) return(x)

  # A one-draw tensor presents its modes directly to the user.
  storage_dims <- dim(unclass(x))
  if (n_draws(x) == 1L && length(args) < length(storage_dims)) {
    args <- c(list(1L), args)
  }

  dims <- storage_dims
  if (length(args) < length(dims)) {
    args <- c(args, rep(list(TRUE), length(dims) - length(args)))
  }

  caller <- parent.frame()
  eval_args <- lapply(args, function(arg) {
    if (is.symbol(arg) && identical(as.character(arg), "")) TRUE
    else eval(arg, caller)
  })

  result <- do.call(`[`, c(list(unclass(x)), eval_args, list(drop = drop)))

  if (!drop && is.array(result) && length(dim(result)) >= 2L) {
    return(.new_tensor_array(result))
  }
  result
}
