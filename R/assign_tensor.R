#' Assign into a tensor
#'
#' Assignment follows ordinary array rules. For a one-draw tensor, the draw
#' axis is implicit and is inserted automatically.
#'
#' @param x A tensor object.
#' @param ... Array indices.
#' @param value A replacement value, including a tensor or ordinary array.
#' @return The modified tensor.
#' @export
#' @method [<- tensor
`[<-.tensor` <- function(x, ..., value) {
  args <- as.list(substitute(list(...)))[-1L]
  if (length(args) == 0L) stop("Indices are required for tensor assignment.", call. = FALSE)
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
  replacement <- if (inherits(value, "tensor")) unclass(value) else value
  out <- do.call(`[<-`, c(list(unclass(x)), eval_args, list(value = replacement)))
  .new_tensor_array(out)
}
