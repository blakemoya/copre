#' Obtain Functionals from a CopRe Result
#'
#' @param obj A \code{copre_result} object.
#' @param f A list of functions.
#' @param mean A logical value indicating whether or not to obtain the
#'  functional from the pointwise mean of the sampled distributions or from each individually.
#'
#' @return The integral over the \code{copre_result} grid of the functions in
#'  the list multiplied by the density of each sample distribution in
#'  \code{obj}.
#' @export
functional <- function(obj, f, mean = FALSE) {
  if (obj$func != 'density') {
    if ('copre_result' %in% class(obj)) {
      obj <- grideval(obj, func = 'density')
    } else {
      stop('`obj` is not a density evaluation, regridify the `mdpolya_result`.')
    }
  }
  if (mean) {
    obj <- matrix(apply(obj, 1, mean), ncol = 1)
  }
  out <- apply(obj, 1, function(y)
      pracma::trapz(obj$grid, f(obj$grid) * y))
  return(out)
}
