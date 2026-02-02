#' tt
#'
#' Computes the tensor-train decomposition of a tensor using the SVD.
#'
#' @param A An array containing the tensor to be decomposed.
#' @param ranks A list of ranks used for the approximation.
#' If provided, accuracy will be ignored.
#' @param accuracy Determines the accuracy of the approximation.
#' If ranks are not provided, will compute them based on \eqn{\delta}-truncation
#'
#' @return A list of the tensor-train cores used to reconstruct the array.
#' @details See I. V. Oseledets, (2011). Tensor-Train Decomposition.
#' SIAM J. SCI. COMPUT. for more details.
#'
#' @examples
#' A <- array(1:24, dim = c(2, 3, 4))
#' res <- tt_svd(A, ranks = c(2, 2, 2))
#'
#' @seealso [tt_reconstruct()] to reconstruct the tensor from the decomposed parts.
#'
#' @export
tt <- function(A, ranks, epsilon, method = "svd") {
  if(missing(ranks) & missing(epsilon)) stop("Provide ranks or epsilon.")

  dims <- dim(A)
  order <- length(dims)

  list_cores <- vector(mode = "list", length = order)
  C <- A

  if(!missing(ranks)) {
    ranks <- c(1, ranks)

    if(length(ranks) != (order+1)) stop("Ranks should be the same length as the
                                        order of A subtracted by 1.")
    for(k in 1:(order-1)) {
      C <- matrix(as.vector(C), nrow = ranks[k] * dims[k])

      C_svd <- svd(C, nu = ranks[k + 1], nv = ranks[k + 1])

      U <- C_svd$u
      d <- C_svd$d
      V <- C_svd$v

      list_cores[[k]] <- drop(array(U, dim = c(ranks[k], dims[k], ranks[k+1])))

      C <- diag(d, nrow = ranks[k+1], ncol = ranks[k+1]) %*% t(V)
    }
  }

  else if(!missing(epsilon)) {
    A_fnorm <- frob_norm(A)

    delta <- (epsilon)/(sqrt(order - 1)) * A_fnorm

    r_delta_past <- 1

    for(k in 1:(order-1)) {
      C <- matrix(as.vector(C), nrow = r_delta_past * dims[k])

      C_svd <- svd(C)

      s <- C_svd$d

      total <- sum(s^2)
      cum_kept <- cumsum(s^2)

      r_delta <- min(which(total - cum_kept <= delta^2)[1], nrow(s))

      U <- C_svd$u[, 1:r_delta, drop = FALSE]
      d <- C_svd$d[1:r_delta]
      V <- C_svd$v[, 1:r_delta, drop = FALSE]

      list_cores[[k]] <- drop(array(U, dim = c(r_delta_past, dims[k], r_delta)))

      C <- diag(d, nrow = r_delta, ncol = r_delta) %*% t(V)

      r_delta_past <- r_delta
    }
  }

  list_cores[[order]] <- C

  list_cores
}

