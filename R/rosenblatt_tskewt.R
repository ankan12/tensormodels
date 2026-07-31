#' Rosenblatt transform for tensor-variate skewed-t observations
#'
#' Convenience wrapper around [tensor_rosenblatt()] for tensor-variate
#' skewed-t observations.
#'
#' @param draws A `tensor` object containing one or more observations.
#' @param mu A tensor/array containing the location parameter.
#' @param skew A tensor/array containing the skew parameter.
#' @param sigmas A list containing one covariance matrix for each tensor mode.
#' @param nu The positive skewed-t degrees-of-freedom parameter.
#' @param rel.tol Relative tolerance passed to [stats::integrate()].
#' @param abs.tol Absolute tolerance passed to [stats::integrate()].
#' @param subdivisions Maximum number of subintervals used by
#'   [stats::integrate()].
#' @param show_progress Whether to report progress after each tensor draw.
#'
#' @return A `tensor` object with the same draw shape and draw count as
#'   `draws`. Under the model with known parameters, its entries are
#'   independent Uniform(0, 1) variables in R's column-major tensor order.
#'
#' @details
#' This wrapper constructs the skewed-t parameter list and calls
#' [tensor_rosenblatt()]. See that function for the conditional GIG
#' quadrature and fitted-parameter caveat.
#'
#' @examples
#' set.seed(1)
#' draws <- rtskewt(
#'   n = 5,
#'   mu = array(0, dim = c(2, 2)),
#'   skew = array(0.2, dim = c(2, 2)),
#'   sigmas = list(diag(2), diag(2)),
#'   nu = 6
#' )
#' uniforms <- rosenblatt_tskewt(
#'   draws,
#'   mu = array(0, dim = c(2, 2)),
#'   skew = array(0.2, dim = c(2, 2)),
#'   sigmas = list(diag(2), diag(2)),
#'   nu = 6
#' )
#' uniforms
#' @export
rosenblatt_tskewt <- function(draws,
                              mu,
                              skew,
                              sigmas,
                              nu,
                              rel.tol = 1e-7,
                              abs.tol = 1e-9,
                              subdivisions = 100L,
                              show_progress = FALSE) {
  tensor_rosenblatt(
    draws = draws,
    model = "skewt",
    parameters = list(
      mu = mu,
      skew = skew,
      sigmas = sigmas,
      nu = nu
    ),
    rel.tol = rel.tol,
    abs.tol = abs.tol,
    subdivisions = subdivisions,
    show_progress = show_progress
  )
}

.tskewt_rosenblatt_conditional_cdf <- function(value,
                                                skew,
                                                lambda,
                                                chi,
                                                psi,
                                                rel.tol,
                                                abs.tol,
                                                subdivisions) {
  .tensor_rosenblatt_gig_cdf(
    value = value,
    skew = skew,
    lambda = lambda,
    chi = chi,
    psi = psi,
    rel.tol = rel.tol,
    abs.tol = abs.tol,
    subdivisions = subdivisions
  )
}
