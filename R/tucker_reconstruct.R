#' tucker_reconstruct
#'
#' Computes the approximated tensor from a Tucker decomposition
#'
#' @param list_cores A list of cores generated from [tucker()]
#'
#' @return The reconstructed tensor generated from the Tucker decomposition.
#' @details See T. Kolda, B. Bader, "Tensor decomposition and applications".
#' SIAM Applied Mathematics and Applications 2009. for more details.
#'
#' @examples
#' A <- array(1:24, dim = c(2, 3, 4))
#' res <- tucker(A, ranks = c(2, 2, 2))
#' tucker_reconstruct(res)
#'
#' @seealso [tucker()] to create the list of cores from the Tucker decomposition.
#'
#' @export
tucker_reconstruct <- function(tucker_cores) {
  G_core <- tucker_cores$G

  mats <- tucker_cores$mats

  order <- length(mats)

  reconstruct <- n_prod(G_core, 1, mats[[1]])

  for(i in 2:order) {
    reconstruct <- n_prod(reconstruct, i, mats[[i]])
  }
  reconstruct
}
