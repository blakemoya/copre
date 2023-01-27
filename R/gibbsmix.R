#' Marginal Gibbs-type Mixture Model Sampler
#'
#' @description A function that samples marginal mixture densities via a
#'  marginal Gibbs sampler.
#'
#' @param data The data from which to sample predictive distributions.
#' @param k The number of predictive samples to draw.
#' @param b_msr A `base_measure` object.
#' @param s_msr A `seq_measure` object.
#' @param burn The number of initial sampling iterations to discard, will be
#'  truncated if a non-integer.
#' @param thin The number of sampling iterations to discard between records,
#'  will be truncated if a non-integer.
#'
#' @return A `seqre_result` object.
#' @seealso [seqre()], [seq_measure()], [base_measure()]
#' @export
gibbsmix <- function(data, k, b_msr, s_msr, burn = 1000, thin = 150) {
  if (!('base_measure' %in% class(b_msr))) {
    stop('`b_msr` must be a `base_measure` object.')
  }
  if (!('seq_measure' %in% class(s_msr))) {
    stop('`s_msr` must be a `seq_measure` object.')
  }
  z0 <- s_msr$rnext(length(data))
  res_marg <- marg_cpp(data, k, z0 - 1,
                       b_msr$idx, b_msr$pars, b_msr$hpars,
                       s_msr$idx, s_msr$pars, s_msr$hpars,
                       burn, thin)

  out <- list()
  out$theta <- res_marg[[1]]
  out$args <- list(data = data, k = k, s_msr = s_msr, b_msr = b_msr,
                   burn = burn, thin = thin)
  phi <- apply(res_marg[[1]], 1, function(t) {
    z <- t[, 1]
    z_uq <- sort(unique(z))
    z_w <- table(z) / length(z)
    mat <- matrix(c(z_w, t[z_uq, 2:(b_msr$dim + 1)]),
                  nrow = length(z_uq), ncol = b_msr$dim + 1)
    return(mat)
  }, simplify = FALSE)
  out$phi <- phi
  class(out) <- 'seqre_result'
  return(out)
}
