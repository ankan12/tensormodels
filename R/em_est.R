#' em_est
#'
#' Estimates the mean array and covariance matrices from an array of tensor
#' variate normal draws using an expectation-maximization algorithm.
#'
#' @param draws An array containing the draws, where the first mode represents
#'   each draw.
#' @param max_iter A max number of iterations to try to get covariance matrices
#'   that converge.
#' @param tol A tolerance level to define the convergence of matrices.
#' @param method Specify the distribution of draws. Must be one of skewt,
#'    vargamma, invgauss, or gen hyper
#' @return A list containing the estimated mean array and the list of covariance
#'   matrices.
#'
#' @export
em_est <- function(draws, max_iter = 1000, tol = 1e-2, quiet = TRUE, model) {
  if (!model %in% c("skewt", "vargamma", "invgauss", "genhyper")) {
    stop("Not a valid model. Must be skewt, vargamma, invgauss, or genhyper")
  }

  if(model == "skewt") em_est_skewt(draws, max_iter, tol, quiet)
  else if(model == "vargamma") em_est_vargamma(draws, max_iter, tol, quiet)
  else if(model == "invgauss") em_est_invgauss(draws, max_iter, tol, quiet)
  else if(model == "genhyper") em_est_genhyper(draws, max_iter, tol, quiet)
}
