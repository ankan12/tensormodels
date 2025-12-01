#' mean_rel_change
#'
#' Computes the mean relative difference of Frobenius norms between parameters.
#' @return A value representing the mean Frobenius difference between all parameters to check for convergence.
#'
#' @noRd
mean_rel_change <- function(mu_new, mu_old, skew_new, skew_old, sig_new, sig_old, nu_new, nu_old) {
  mu_frob <- frob_norm_diff(mu_new, mu_old)
  skew_frob <- frob_norm_diff(skew_new, skew_old)

  sig_frobs <- rep(0, length(sig_new))

  for(k in 1:length(sig_new)) {
    sig_frobs[k] <- frob_norm_diff(sig_new[[k]], sig_old[[k]])
  }

  nu_frob <- frob_norm_diff(nu_new, nu_old)

  mean(c(mu_frob, skew_frob, sig_frobs, nu_frob))
}
