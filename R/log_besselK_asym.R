.log_besselK_asym <- function(x, nu) {
  Bessel::besselK.nuAsym(
    x = x,
    nu = abs(nu),
    k.max = 5,
    expon.scaled = TRUE,
    log = TRUE
  )
}

.besselK_asym_ratio <- function(x, numerator_nu, denominator_nu) {
  exp(.log_besselK_asym(x, numerator_nu) -
        .log_besselK_asym(x, denominator_nu))
}

.dlog_besselK_asym_dnu <- function(x, nu, eps = 1e-5) {
  (.log_besselK_asym(x, nu + eps) -
     .log_besselK_asym(x, nu - eps)) / (2 * eps)
}
