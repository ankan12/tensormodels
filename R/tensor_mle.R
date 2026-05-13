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
#' @param model Specify the distribution of draws. Must be one of normal, skewt,
#'    vargamma, invgauss, or gen hyper
#' @param penalize_nu Logical. For the skewt model, should a shifted gamma
#'   prior penalty be applied to the degrees of freedom?
#' @param nu_prior_shape,nu_prior_rate Shape and rate for the shifted gamma
#'   prior \code{nu - 4 ~ Gamma(shape, rate)} used when \code{penalize_nu}
#'   is \code{TRUE}.
#' @return A list containing the estimated mean array and the list of covariance
#'   matrices.
#'
#' @export
tensor_mle <- function(draws, max_iter = 1000, tol = 1e-2, quiet = TRUE,
                       model, restrict = NULL,
                       penalize_nu = FALSE, nu_prior_shape = 2,
                       nu_prior_rate = 0.05) {
  if (!model %in% c("normal", "skewt", "vargamma", "invgauss", "genhyper")) {
    stop("Not a valid model. Must be normal, skewt, vargammax, invgauss, or genhyper")
  } # check there is a real model

  if (!is.list(draws) || length(draws) == 0) {
    stop("Draws must be a non-empty list.")
  } # check the draws are in a list

  n <- length(draws)

  dim_first <- dim(draws[[1]])

  for(i in 2:n) { # check all draws are same dimension
    curr_dim <- dim(draws[[i]])

    if(!all.equal(dim_first, curr_dim)) {
        stop(sprintf(
        "All draws must have the same dimensions. Draw 1 has dimensions %s, but draw %d has dimensions %s.",
        paste(dim_first, collapse = " x "), i,
        paste(curr_dim, collapse = " x "))
      )
    }
  }

  # call the correct MLE function based on model
  if(model == "normal") tensor_mle_normal(draws, max_iter, tol, quiet, restrict)
  else if(model == "skewt") {
    tensor_mle_skewt(draws, max_iter, tol, quiet, restrict,
                     penalize_nu = penalize_nu,
                     nu_prior_shape = nu_prior_shape,
                     nu_prior_rate = nu_prior_rate)
  }
  else if(model == "vargamma") tensor_mle_vargamma(draws, max_iter, tol, quiet, restrict)
  else if(model == "invgauss") tensor_mle_invgauss(draws, max_iter, tol, quiet, restrict)
  else if(model == "genhyper") tensor_mle_genhyper(draws, max_iter, tol, quiet, restrict)
}
