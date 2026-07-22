#' Construct tensor-valued observations
#'
#' `tensor()` stores one or more identically shaped tensor observations. By
#' default an input array is interpreted as one draw; `obs` identifies an
#' existing observation dimension when an array contains multiple draws.
#'
#' @param data An array, or the data used to construct one or more draws when
#'   `dim` is supplied.
#' @param dim Optional dimensions of each draw, following [array()].
#' @param obs `NULL` for one draw, or a numeric/character observation axis.
#' @param n Number of draws to construct when `dim` is supplied. Defaults to 1.
#' @return A `tensor` object with observations on its first dimension.
#' @export
tensor <- function(data, dim = NULL, obs = NULL, n = 1L, dimnames = NULL) {
  if (inherits(data, "tensor")) {
    if (!is.null(dim) || !is.null(obs) || n != 1L || !is.null(dimnames)) {
      stop("Additional constructor arguments must be NULL/default when `data` is already a tensor.", call. = FALSE)
    }
    return(data)
  }

  if (length(n) != 1L || is.na(n) || !is.numeric(n) || !is.finite(n) ||
      n != floor(n) || n < 1L) {
    stop("`n` must be one positive integer.", call. = FALSE)
  }
  n <- as.integer(n)

  if (!is.null(dim)) {
    if (!is.numeric(dim) || length(dim) == 0L || anyNA(dim) ||
        any(!is.finite(dim)) || any(dim < 1) || any(dim != floor(dim))) {
      stop("`dim` must be a vector of positive integers.", call. = FALSE)
    }
    dim <- as.integer(dim)
    expected <- prod(dim) * n
    if (length(data) != expected) {
      warning(
        sprintf("Data length %d does not match requested tensor size %d; values will be recycled.",
                length(data), expected),
        call. = FALSE
      )
    }
    if (!is.null(obs)) stop("`obs` cannot be supplied with `dim`.", call. = FALSE)

    # Construct in ordinary array order (one complete draw at a time), then
    # move the draw axis to the front for the package's internal layout.
    if (is.null(dimnames)) {
      raw <- array(data, dim = c(dim, n))
    } else {
      raw <- array(data, dim = c(dim, n),
                   dimnames = c(dimnames, list(NULL)))
    }
    values <- aperm(raw, c(length(dim) + 1L, seq_along(dim)))
    if (n == 1L && !is.null(dimnames)) {
      dimnames(values) <- c(list(NULL), dimnames)
    }
    return(.new_tensor_array(values))
  }

  if (!is.array(data)) stop("`data` must be an array when `dim` is NULL.", call. = FALSE)
  if (n != 1L) stop("`n` is only available when `dim` is supplied.", call. = FALSE)

  if (!is.null(obs)) {
    obs <- vctrs::vec_as_location2(obs, n = length(dim(data)),
                                   names = names(dim(data)), arg = "obs")
    if (obs != 1L) data <- aperm(data, c(obs, setdiff(seq_along(dim(data)), obs)))
    return(.new_tensor_array(unclass(data)))
  }

  if (is.null(dimnames(data))) {
    values <- array(unclass(data), dim = c(1L, dim(data)))
  } else {
    values <- array(unclass(data), dim = c(1L, dim(data)),
                    dimnames = c(list(NULL), dimnames(data)))
  }
  .new_tensor_array(values)
}

.new_tensor_array <- function(x) {
  if (!is.array(x) || length(dim(x)) < 2L) {
    stop("Internal tensor storage must have a draw axis and tensor modes.", call. = FALSE)
  }
  structure(unclass(x), class = c("tensor", "array"))
}

.tensor_slice_array <- function(x, i) {
  do.call(`[`, c(list(unclass(x), i), rep(list(TRUE), length(draw_shape(x))),
                 list(drop = FALSE)))
}

.tensor_broadcast_array <- function(x, n) {
  if (n_draws(x) == n) return(unclass(x))
  if (n_draws(x) != 1L) {
    stop("Only a one-draw tensor can be broadcast to another draw count.", call. = FALSE)
  }
  .tensor_slice_array(x, rep.int(1L, n))
}

.tensor_single_draw_array <- function(x, arg = "x") {
  if (!inherits(x, "tensor")) return(x)
  if (n_draws(x) != 1L) stop(sprintf("`%s` must contain exactly one draw.", arg), call. = FALSE)
  out <- unclass(x)
  dn <- dimnames(out)[-1L]
  dim(out) <- draw_shape(x)
  if (!is.null(dn)) dimnames(out) <- dn
  out
}

#' Arithmetic and comparison operations for tensors
#'
#' Operations are applied to every draw. A one-draw tensor is broadcast across
#' a tensor with more draws.
#' @param e1,e2 Tensor or numeric operands.
#' @return Arithmetic operations return a tensor; comparisons return an array.
#' @export
Ops.tensor <- function(e1, e2) {
  arithmetic <- c("+", "-", "*", "/", "^", "%%", "%/%")
  e1_tensor <- inherits(e1, "tensor")
  if (missing(e2)) {
    out <- do.call(.Generic, list(unclass(e1)))
    return(if (.Generic %in% arithmetic) .new_tensor_array(out) else out)
  }
  e2_tensor <- inherits(e2, "tensor")

  if (e1_tensor && e2_tensor) {
    if (!identical(unname(draw_shape(e1)), unname(draw_shape(e2)))) {
      stop("Tensor draw shapes must agree.", call. = FALSE)
    }
    n1 <- n_draws(e1); n2 <- n_draws(e2)
    if (n1 != n2 && n1 != 1L && n2 != 1L) {
      stop("Tensor draw counts must agree or one must equal 1.", call. = FALSE)
    }
    n <- max(n1, n2)
    lhs <- .tensor_broadcast_array(e1, n)
    rhs <- .tensor_broadcast_array(e2, n)
  } else {
    tx <- if (e1_tensor) e1 else e2
    other <- if (e1_tensor) e2 else e1
    lhs <- if (e1_tensor) unclass(e1) else other
    rhs <- if (e2_tensor) unclass(e2) else other
    if (length(other) != 1L && length(other) != prod(draw_shape(tx))) {
      stop("A non-tensor operand must have length 1 or one value per tensor entry.", call. = FALSE)
    }
    if (length(other) > 1L) {
      shaped <- array(other, dim = c(1L, draw_shape(tx)))
      shaped <- .tensor_slice_array(.new_tensor_array(shaped), rep.int(1L, n_draws(tx)))
      if (e1_tensor) rhs <- shaped else lhs <- shaped
    }
  }
  out <- do.call(.Generic, list(lhs, rhs))
  if (.Generic %in% arithmetic) .new_tensor_array(out) else out
}
