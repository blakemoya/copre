#' Obtain Moments from a CopRe or SeqRe Result
#'
#' @param obj A `copre_result` or `seqre_result` object.
#' @param mom A numeric scalar indicating the moment to calculate.
#' @param cntrl A logical value indicating whether the moment should be central
#'   or not. Defaults to `TRUE`.
#' @param grd A numeric vector of grid values on which the density function
#'   samples in `obj` should be calculated for trapezoidal integration.
#'
#' @return A vector of moment values for each sampled distribution in `obj`.
#' @export
moment <- function(obj, mom, cntrl = TRUE, grd = NULL) {
  UseMethod('moment')
}

#' @describeIn moment Moment calculation method for `seqre_result` objects.
#' @export
moment.seqre_result <- function(obj, mom, cntrl = TRUE, grd = NULL) {
  moments(grideval(obj, grd = grd, func = 'density'), mom = mom, cntrl = cntrl)
}

#' @describeIn moment Moment calculation method for `grideval_result` objects.
#' @export
moment.grideval_result <- function(obj, mom, cntrl = TRUE, grd = NULL) {
  if (obj$func != 'density') {
    if ('copre_result' %in% class(obj)) {
      obj <- grideval(obj, func = 'density')
    } else {
      stop('`obj` is not a density evaluation, regridify the `mdpolya_result`.')
    }
  }
  x <- obj$grid
  x_cntr <- x
  if (cntrl & (mom != 1)) {
    x_cntr <- obj$grid - apply(obj, 1, function(f)
      pracma::trapz(obj$grid, obj$grid * f))
  }
  return(apply(obj, 1, function(f) pracma::trapz(x, x_cntr ^ mom * f)))
}
