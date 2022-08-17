#' @param obj An \code{mdpolya_result} object.
#' @param i A numeric vector of sample indices.
#'
#' @describeIn mdp Subset method for \code{mdpolya_result} objects
#' @export
`[[.mdpolya_result` <- function(obj, i) {
  if (!is.null(obj$theta)) {
    obj$theta <- obj$theta[i, , , drop = FALSE]
  }
  obj$eta <- obj$eta[i, ]
  obj$phi <- obj$phi[i]
  return(obj)
}

#' @param obj A \code{grideval_result} object.
#' @param name The name of the attribute to access (i.e. \code{func},
#'   \code{grid}, or \code{args}).
#'
#' @describeIn grideval Attribute access method for \code{grideval_result}
#'   objects
#' @export
`$.grideval_result` <- function(obj, name) {
  attr(obj, name)
}

#' @param obj A \code{grideval_result} object.
#' @param i A numeric vector of sample indices.
#'
#' @describeIn grideval Subset method for \code{grideval_result} objects
#' @export
`[[.grideval_result` <- function(obj, i) {
  out <- matrix(obj, dim(obj))[i, , drop = FALSE]
  class(out) <- class(obj)
  attr(out, 'func') <- attr(obj, 'func')
  attr(out, 'grid') <- attr(obj, 'grid')
  attr(out, 'args') <- attr(obj, 'args')
  attr(out, 'args')$k <- length(i)
  return(out)
}
