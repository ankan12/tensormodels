#' tensor_mle
#'
#' Estimates the mean array and covariance matrices from an array of tensor
#' variate draws. If tensor variate normal, uses the flop-flop algorithm. Otherwise,
#' uses the expectation-maximization algorithm.
#'
#' @param draws An array containing the draws, where the first mode represents
#'   each draw.
#' @param max_iter A max number of iterations to try to get covariance matrices
#'   that converge.
#' @param tol A tolerance level to define the convergence of matrices.
#' @param method Specify the distribution of draws. Must be one of normal, skewt,
#'    vargamma, invgauss, or gen hyper
#' @return A list containing the estimated mean array and the list of covariance
#'   matrices.
#'
#' @export
tensor_mle <- function(draws, max_iter = 1000, tol = 1e-2, quiet = TRUE, model) {
  if (!model %in% c("normal", "skewt", "vargamma", "invgauss", "genhyper")) {
    stop("Not a valid model. Must be normal, skewt, vargamma, invgauss, or genhyper")
  }

  if(model == "normal") tensor_mle_normal(draws, max_iter, tol, quiet)
  else if(model == "skewt") tensor_mle_skewt(draws, max_iter, tol, quiet)
  else if(model == "vargamma") tensor_mle_vargamma(draws, max_iter, tol, quiet)
  else if(model == "invgauss") tensor_mle_invgauss(draws, max_iter, tol, quiet)
  else if(model == "genhyper") tensor_mle_genhyper(draws, max_iter, tol, quiet)
}
