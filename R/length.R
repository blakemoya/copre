#' Length
#'
#' @param x A `grideval_result` object.
#'
#' @return The number of samples `k` in `obj`.
#' @export
length.grideval_result <- function(x) {
  return(nrow(x))
}
