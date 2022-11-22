#' Sequence Resampling
#'
#' @description A function that samples predictive distributions for univariate
#'  continuous data using exchangeable predictive extension.
#'
#' @param obj A \code{seqre_result} object, usually output from
#'  \code{gibbsmix()}.
#' @param inc A positive integer increment value for the number of predictive
#'  samples to take each convergence check.
#' @param eps An error value which determines the convergence approximation.
#' @param max_it A positive integer maximum number of iterations before halting.
#'
#' @return A \code{seqre_result} object, or a list of two \code{seqre_result}
#'  objects if \code{keep_marg} is \code{TRUE}.
#' @seealso [gibbsmix()]
#' @export
seqre <- function(obj, inc = 1000, eps = 0.001, max_it = 100) {
  data <- obj$args$data
  b_msr <- obj$args$b_msr
  s_msr <- obj$args$s_msr
  if (!('base_measure' %in% class(b_msr))) {
    stop('`b_msr` must be a `base_measure` object.')
  }
  if (!('seq_measure' %in% class(s_msr))) {
    stop('`s_msr` must be a `seq_measure` object.')
  }
  phi <- seqre_cpp(obj$phi, length(data),
                   b_msr$idx, b_msr$pars, b_msr$hpars,
                   s_msr$idx, s_msr$pars, s_msr$hpars,
                   eps, inc, max_it)
  obj$phi <- phi
  obj$args$inc <- inc
  obj$args$eps <- eps
  obj$args$max_it <- max_it
  return(obj)
}
