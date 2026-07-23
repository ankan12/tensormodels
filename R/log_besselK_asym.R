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
  ratio <- exp(.log_besselK_asym(x, numerator_nu) -
                 .log_besselK_asym(x, denominator_nu))

  # The asymptotic implementation is undefined at some small orders, notably
  # K_0(x) in the scalar inverse-Gaussian E-step. Exponential scaling cancels
  # in the ratio and avoids underflow in the direct fallback.
  bad <- !is.finite(ratio)
  if (any(bad)) {
    ratio[bad] <- besselK(
      x[bad], nu = abs(numerator_nu), expon.scaled = TRUE
    ) / besselK(
      x[bad], nu = abs(denominator_nu), expon.scaled = TRUE
    )
  }

  ratio
}

.dlog_besselK_asym_dnu <- function(x, nu, eps = 1e-5) {
  (.log_besselK_asym(x, nu + eps) -
     .log_besselK_asym(x, nu - eps)) / (2 * eps)
}
