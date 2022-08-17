#' Grid evaluation of \code{copre_result} and \code{mdpolya_result} objects
#'
#' @param obj A \code{copre_result} or \code{mdpolya_result} object.
#' @param grd For \code{mdpolya_result} objects, a numeric vector of \code{m}
#'  grid points.
#' @param func Either 'distribution', 'density', or 'gradient'.
#' @param nthreads The number of parallel threads to launch with OpenMP.
#'
#' @return A \code{grideval_result} object, which is a matrix with dimension
#'  \code{[k, m]} of evaluated sample functions, with the following attributes:
#'  * \code{func}: The evaluated function.
#'  * \code{grid}: The grid points on which each of the \code{k} rows was
#'   evaluated.
#'  * \code{args}: A copy of the \code{args} entry from \code{obj}.
#' @export
grideval <- function(obj, grd = NULL, func = 'density', nthreads = 1) {
  UseMethod('grideval')
}

#' @describeIn grideval Grid evaluation method for \code{copre_result} objects.
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

#' @describeIn grideval Grid evaluation method for \code{mdpolya_result}
#'  objects.
#' @export
grideval.mdpolya_result <- function(obj, grd = NULL, func = 'density',
                                    nthreads = 1) {
  if (func == 'distribution') {
    f <- 0
  } else if (func == 'density') {
    f <- 1
  } else if (func == 'gradient') {
    f <- 2
  } else {
    stop(paste('Unrecognized `func`. Choose either \'distribution\',',
               '\'density\' or \'gradient\'.'))
  }
  if (is.null(grd)) {
    r_x <- range(obj$args$data)
    rr_x <- diff(r_x) / 10
    grd <- seq(r_x[1] - rr_x, r_x[2] + rr_x, length = 1000)
  }
  out <- evalmdpolya_cpp(obj$phi, grd, f, nthreads)
  attr(out, 'func') <- func
  attr(out, 'grid') <- grd
  attr(out, 'args') <- obj$args
  class(out) <- 'grideval_result'
  return(out)
}
