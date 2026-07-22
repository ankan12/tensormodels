#' Structure of a tensor
#' @param object A tensor object.
#' @param ... Unused.
#' @return A compact description.
#' @export
str.tensor <- function(object, ...) {
  sprintf("%d tensor %s of shape %s <%s>", n_draws(object),
          if (n_draws(object) == 1L) "draw" else "draws",
          paste(draw_shape(object), collapse = " x "), typeof(unclass(object)))
}
