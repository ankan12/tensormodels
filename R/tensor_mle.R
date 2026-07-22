#' tensor_mle
#'
#' Estimates the mean array and covariance matrices from tensor-valued
#' observations. If tensor variate normal, uses the flop-flop algorithm. Otherwise,
#' uses the expectation-maximization algorithm.
#'
#' @param draws A `tensor` object containing IID observations.
#' @param max_iter A max number of iterations to try to get covariance matrices
#'   that converge.
#' @param tol A tolerance level to define the convergence of matrices.
#' @param model Specify the distribution of draws. Must be one of normal, skewt,
#'    vargamma, invgauss, or genhyper
#' @return A list containing the estimated mean array and the list of covariance
#'   matrices.
#'
#' @export
tensor_mle <- function(draws, max_iter = 1000, tol = 1e-4, quiet = TRUE,
                       model, restrict = NULL) {
  if (!inherits(draws, "tensor")) {
    stop(
      "`draws` must be a `tensor` object. Use `tensor()` first.",
      call. = FALSE
    )
  }

  if (!model %in% c("normal", "skewt", "vargamma", "invgauss", "genhyper")) {
    stop("Not a valid model. Must be normal, skewt, vargamma, invgauss, or genhyper")
  } # check there is a real model

  n <- n_draws(draws)

  if (n == 0L) {
    stop("`draws` must contain at least one observation.", call. = FALSE)
  }

  # call the correct MLE function based on model
  if(model == "normal") tensor_mle_normal(draws, max_iter, tol, quiet, restrict)
  else if(model == "skewt") tensor_mle_skewt(draws, max_iter, tol, quiet, restrict)
  else if(model == "vargamma") tensor_mle_vargamma(draws, max_iter, tol, quiet, restrict)
  else if(model == "invgauss") tensor_mle_invgauss(draws, max_iter, tol, quiet, restrict)
  else if(model == "genhyper") tensor_mle_genhyper(draws, max_iter, tol, quiet, restrict)
}
