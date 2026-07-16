#' Construct a tensor
#'
#' `tensor()` wraps a base array as a tensor, or constructs an array from raw
#' data when `dim` is supplied. Every array dimension represents a tensor mode.
#'
#' @param x An array, or the data used to construct one.
#' @param dim An optional integer vector giving the dimensions of the tensor.
#'   If `NULL`, `x` must already be an array.
#' @param dimnames An optional list of dimension names used when `dim` is
#'   supplied.
#'
#' @return An object of class `tensor` backed by `x`.
#'
#' @examples
#' A <- tensor(array(1:24, dim = c(2, 3, 4)))
#' B <- tensor(1:24, dim = c(2, 3, 4))
#' A + A
#' 2 * A
#'
#' @export
tensor <- function(x, dim = NULL, dimnames = NULL) {
  if (is.null(dim)) {
    if (!is.array(x)) {
      stop(
        "Supply `dim` when `x` is not already an array.",
        call. = FALSE
      )
    }

    if (!is.null(dimnames)) {
      dimnames(x) <- dimnames
    }
  } else {
    if (!is.numeric(dim) || length(dim) == 0L || anyNA(dim) ||
        any(!is.finite(dim)) || any(dim < 1) || any(dim != floor(dim)) ||
        any(dim > .Machine$integer.max)) {
      stop("`dim` must be a vector of positive integers.", call. = FALSE)
    }

    dim <- as.integer(dim)

    if (length(x) != prod(dim)) {
      stop(
        "The length of `x` must equal the product of `dim`.",
        call. = FALSE
      )
    }

    x <- array(x, dim = dim, dimnames = dimnames)
  }

  structure(
    unclass(x),
    class = c("tensor", "array")
  )
}

#' Arithmetic and comparison operations for tensors
#'
#' Tensor arithmetic is elementwise and preserves the tensor class. Operations
#' between two tensors require identical shapes. Comparisons and logical
#' operations return ordinary arrays.
#'
#' @param e1,e2 Tensor or numeric operands.
#'
#' @return Arithmetic operations return a tensor. Comparison and logical
#'   operations return a base array.
#' @export
Ops.tensor <- function(e1, e2) {
  arithmetic <- c("+", "-", "*", "/", "^", "%%", "%/%")
  e1_is_tensor <- inherits(e1, "tensor")

  if (missing(e2)) {
    result <- do.call(.Generic, list(unclass(e1)))

    if (.Generic %in% arithmetic && is.array(result)) {
      return(tensor(result))
    }

    return(result)
  }

  e2_is_tensor <- inherits(e2, "tensor")

  if (
    e1_is_tensor && e2_is_tensor &&
      !identical(dim(e1), dim(e2))
  ) {
    stop(
      sprintf(
        "Tensor shapes must agree: %s and %s.",
        paste(dim(e1), collapse = " x "),
        paste(dim(e2), collapse = " x ")
      ),
      call. = FALSE
    )
  }

  lhs <- if (e1_is_tensor) unclass(e1) else e1
  rhs <- if (e2_is_tensor) unclass(e2) else e2
  result <- do.call(.Generic, list(lhs, rhs))

  if (.Generic %in% arithmetic && is.array(result)) {
    return(tensor(result))
  }

  result
}
