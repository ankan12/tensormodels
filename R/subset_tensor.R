#' Subset a tensor.
#'
#' @param A A tensor object (S3 wrapper around an array).
#' @param ... Indices.
#' @param drop Logical; passed to base `[` (defaults to FALSE).
#'
#' @return A tensor that was subset.
#'
#' @examples
#' A <- tensor(1:24, dim = c(2, 3, 4))
#' A[1, ]
#' @export
#' @method [ tensor
`[.tensor` <- function(x, ..., drop = TRUE) {
  dims <- dim(x) #dimensions of input

    args <- as.list(substitute(list(...)))[-1] #args provided

  # fill missing with TRUE
  if (length(args) < length(dims)) {
    args <- c(args, rep(list(TRUE), length(dims) - length(args)))
  }

  # replace empty arguments (e.g., `A[1, ]`) with TRUE
  for (i in seq_along(args)) {
    if (length(args[[i]]) == 0) args[[i]] <- TRUE
  }

  # Evaluate indices safely
  eval_indices <- lapply(args, function(a) {
    if (is.symbol(a) && as.character(a) == "") TRUE else eval(a, parent.frame())
  })

  # subset the underlying array
  result <- do.call(`[`, c(list(unclass(x)), eval_indices, list(drop = drop)))

  # return as tensor class
  if (is.array(result)) class(result) <- "tensor"
  result
}


  #full_index <- c(list(i), lapply(dims-1, function(d) TRUE))
  # Detect which args are actually missing
  # is_missing <- vapply(args, function(a) length(a) == 0, logical(1))
  #
  # # Debug prints (optional)
  # print(args)
  # print(is_missing)
  #
  # # Condition: only the first index is provided, all others missing
  # if (length(args) == 1) { # && all(is_missing[-1])) {
  #   # Evaluate first index safely
  #   i <- eval(args[[1]], parent.frame())
  #
  #
  #   # Construct a full index vector for first dim, all others = TRUE
  #   full_index <- c(list(i), lapply(dims-1, function(d) TRUE))
  #
  #   # Subset safely
  #   return(do.call(`[`, c(list(x$data), full_index, list(drop = drop))))
  #
  #   # Slice along first dimension, take all for remaining dims
  #   # dims <- length(dim(x))
  #   # call_args <- c(list(i), rep(list(quote(expr = )), dims - 1))
  #   # return(do.call(`[`, c(list(x), call_args, list(drop = drop))))
  # }
  #
  # # Otherwise, normal subsetting
  # eval_indices <- lapply(args, function(a) eval(a, parent.frame()))
  # do.call(`[`, c(list(x), eval_indices, list(drop = drop)))



