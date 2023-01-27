#' @param obj An `seqre_result` object.
#' @param i A numeric vector of sample indices.
#'
#' @describeIn seqre Subset method for `seqre_result` objects
#' @export
`[[.seqreresult` <- function(obj, i) {
  if (!is.null(obj$theta)) {
    obj$theta <- obj$theta[i, , , drop = FALSE]
  }
  obj$eta <- obj$eta[i, ]
  obj$phi <- obj$phi[i]
  return(obj)
}

#' @param obj A `grideval_result` object.
#' @param name The name of the attribute to access (i.e. `func`, `grid`, or
#'   `args`).
#'
#' @describeIn grideval Attribute access method for `grideval_result` objects
#' @export
`$.grideval_result` <- function(obj, name) {
  attr(obj, name)
}

#' @param obj A `grideval_result` object.
#' @param i A numeric vector of sample indices.
#'
#' @describeIn grideval Subset method for `grideval_result` objects
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
