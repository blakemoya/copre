#' Grid evaluation of `copre_result` and `seqre_result` objects
#'
#' @param obj A `copre_result` or `seqre_result` object.
#' @param grd For `seqre_result` objects, a numeric vector of `m` grid points.
#' @param func Either 'distribution', 'density', or 'gradient'.
#' @param nthreads The number of parallel threads to launch with OpenMP.
#'
#' @return A `grideval_result` object, which is a matrix with dimension `[k, m]`
#'   of evaluated sample functions, with the following attributes:
#'  * `func`: The evaluated function.
#'  * `grid`: The grid points on which each of the `k` rows was evaluated.
#'  * `args`: A copy of the `args` entry from `obj`.
#' @export
grideval <- function(obj, grd = NULL, func = 'density', nthreads = 1) {
  UseMethod('grideval')
}

#' @describeIn grideval Grid evaluation method for `copre_result` objects.
#' @export
grideval.copre_result <- function(obj, grd = NULL, func = 'density',
                                  nthreads = 1) {
  fdict <- c('distribution' = 0, 'density' = 1, 'gradient' = 2)
  if (!(fdict[func] %in% 0:2)) {
    stop(paste('Unrecognized `func`. Choose either \'distribution\',',
               '\'density\' or \'gradient\'.'))
  }
  if (fdict[obj$func] == fdict[func]) {
    return(obj)
  } else if (fdict[obj$func] < fdict[func]) {
    core <- function(obj) {
      out <- t(apply(obj, 1, pracma::gradient, h1 = obj$grid))
      attr(out, 'func') <- names(fdict)[fdict[obj$func] + 2]
      attr(out, 'grid') <- obj$grid
      attr(out, 'args') <- obj$args
      class(out) <- class(obj)
      return(out)
    }
  } else {
    stop('`grideval` should not be used for integration of a `copre_result`.')
    core <- function(obj) {
      trapz <- function(y, x) {pracma::cumtrapz(x, y)}
      out <- t(apply(obj, 1, trapz, x = obj$grid))
      attr(out, 'func') <- names(fdict)[fdict[obj$func]]
      attr(out, 'grid') <- obj$grid
      attr(out, 'args') <- obj$args
      class(out) <- class(obj)
      return(out)
    }
  }
  out <- obj
  while (fdict[out$func] != fdict[func]) {
    out <- core(out)
  }
  return(out)
}

#' @describeIn grideval Grid evaluation method for `seqre_result`
#'  objects.
#' @export
grideval.seqre_result <- function(obj, grd = NULL, func = 'density',
                                  nthreads = 1) {
  if (is.null(grd)) {
    r_x <- range(obj$args$data)
    rr_x <- diff(r_x) / 10
    grd <- seq(r_x[1] - rr_x, r_x[2] + rr_x, length = 1000)
  }
  out <- obj$args$b_msr$eval(obj$phi, grd, func, nthreads)
  attr(out, 'func') <- func
  attr(out, 'grid') <- grd
  attr(out, 'args') <- obj$args
  class(out) <- 'grideval_result'
  return(out)
}
