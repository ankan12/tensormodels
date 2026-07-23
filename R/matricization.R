#' matricization
#'
#' Compute the matricization of a tensor along the kth mode.
#'
#' @param array An array representing a tensor.
#' @param k A mode to matricize along.
#'
#' @return An array: the result of the nm-mode product between A and B.
#'
#' @examples
#' X <- array(1:24, dim = c(2, 3, 4))
#' matricization(X, 1)
#' matricization(X, 2)
#' matricization(X, 3)
#' @export
matricization <- function(X, k) {
  if (inherits(X, "tensor")) {
    shape <- draw_shape(X)
    if (k < 1 || k > length(shape)) stop("Invalid mode index k.")

    mats <- lapply(seq_len(n_draws(X)), function(i) {
      matricization(.tensor_single_draw_array(pull_draw(X, i)), k)
    })

    out_dim <- dim(mats[[1L]])
    out <- array(unlist(mats, use.names = FALSE),
                 dim = c(out_dim, length(mats)))
    return(tensor(out, obs = 3L))
  }

  # Check tensor structure
  if (is.null(dim(X))) return(X)

  dims <- dim(X)
  N <- length(dims)
  if (k < 1 || k > N) stop("Invalid mode index k.")

  # Compute the permutation: bring mode k to the front
  perm <- c(k, setdiff(seq_len(N), k))
  X_perm <- aperm(X, perm)

  # Reshape to matrix
  I_k <- dims[k]
  J <- prod(dims[-k])
  matrix(X_perm, nrow = I_k, ncol = J)
}
