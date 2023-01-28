#' Mode Extractor
#'
#' @description Extracts the modes from a `copre_result` or `seqre_result`
#'   object.
#'
#' @param obj A `copre_result` or `mdp_result` object.
#' @param mean A logical value indicating whether to extract the modes of the
#'  mean density of each of the individual sampled density.
#' @param grd For a `seqre_result`, a grid on which to evaluate the object.
#' @param idx A logical value indicating whether to also return the index within
#'  `grd` of the discovered modes.
#' @param anti A logical value indicating whether to extract true modes or
#'  anti-modes.
#'
#' @return A matrix of modes values in the support of the `copre_result` density
#' @export
modes <- function(obj, mean = FALSE, grd = NULL, idx = FALSE, anti = FALSE) {
  UseMethod('modes')
}

#' Antiode Extractor
#'
#' @description Extracts the antimodes from a `copre_result` or `seqre_result`
#'  object.
#'
#' @param obj A `copre_result` or `mdp_result` object.
#' @param mean A logical value indicating whether to extract the modes of the
#'  mean density of each of the individual sampled density.
#' @param grd For `mdpolya_result`, a grid on which to evaluate the object.
#' @param idx A logical value indicating whether to also return the index within
#'  `grd` of the discovered modes.
#'
#' @return A matrix of antimodes values in the support of the `copre_result`
#'  density
#' @export
antimodes <- function(obj, mean = FALSE, grd = NULL, idx = FALSE) {
  modes(obj, mean = mean, grd = grd, idx = idx, anti = TRUE)
}

#' @describeIn modes Mode-counting method for `seqre_result` objects.
#' @export
modes.seqre_result <- function(obj, mean = FALSE, grd = NULL, idx = FALSE,
                               anti = FALSE) {
  modes(grideval(obj, grd = grd, func = 'density'), mean = mean, anti = anti)
}

#' @describeIn modes Mode-counting method for `grideval_result` objects.
#' @export
modes.grideval_result <- function(obj, mean = FALSE, grd = NULL, idx = FALSE,
                                  anti = FALSE) {
  if (obj$func != 'density') {
    if ('copre_result' %in% class(obj)) {
      obj <- grideval(obj, func = 'density')
    } else {
      stop('`obj` is not a density evaluation, regridify the `mdpolya_result`.')
    }
  }
  if (anti) {
    obj <- -obj
  }
  if (idx) {
    f <- function(p) list(idx = pracma::findpeaks(p)[, 2],
                          value = obj$grid[pracma::findpeaks(p)[, 2]])
  } else {
    f <- function(p) obj$grid[pracma::findpeaks(p)[, 2]]
  }
  if (mean) {
    return(f(apply(obj, 2, mean)))
  } else {
    return(apply(obj, 1, f, simplify = FALSE))
  }
}

#' @param obj A `copre_result` or `seqre_result` object.
#' @param mean A logical value indicating whether to count the modes of the mean
#'  density of each of the individual sampled density.
#' @param grd For `seqre_result`, a grid on which to evaluate the object.
#' @param anti A logical value indicating whether to extract true modes or
#'  anti-modes (i.e. local minima of th3e density function).
#'
#' @describeIn modes Counts the modes from a `copre_result` or `seqre_result`
#'   object.
#' @export
n_modes <- function(obj, mean = FALSE, grd = NULL, anti = FALSE) {
  if (mean) {
    return(length(modes(obj, mean = TRUE, grd = grd, anti = anti)))
  } else {
    return(sapply(modes(obj, mean = FALSE, grd = grd, anti = anti), length))
  }
}


