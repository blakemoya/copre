#' Sequence Resampling
#'
#' @description A function that samples predictive distributions for univariate
#'  continuous data using exchangeable predictive extension.
#'
#' @param data The data from which to sample predictive distributions.
#' @param k The number of predictive samples to draw.
#' @param b_msr A \code{base_measure} object.
#' @param s_msr A \code{seq_measure} object.
#' @param keep_marg A logical value indicating whether or not to keep the
#'  results from the marginal MCMC sampler.
#' @param burn The number of initial sampling iterations to discard, will be
#'  truncated if a non-integer.
#' @param thin The number of sampling iterations to discard between records,
#'  will be truncated if a non-integer.
#' @param inc A positive integer increment value for the number of predictive
#'  samples to take each convergence check.
#' @param eps An error value which determines the convergence approximation.
#' @param max_it A positive integer maximum number of iterations before halting.
#'
#' @return A \code{seqre_result} object, or a list of two \code{seqre_result}
#'  objects if \code{keep_marg} is \code{TRUE}.
#' @seealso [seq_measure()], [base_measure()]
#' @export
seqre <- function(data, k, b_msr, s_msr, keep_marg = FALSE,
                  burn = 1000, thin = 150, inc = 1000, eps = 0.01,
                  max_it = 100) {
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
                   keep_marg = keep_marg, burn = burn, thin = thin)
  phi0 <- apply(res_marg[[1]], 1, function(t) {
    z <- t[, 1]
    z_uq <- sort(unique(z))
    z_w <- table(z) / length(z)
    mat <- matrix(c(z_w, t[z_uq, 2:(b_msr$dim + 1)]),
                  nrow = length(z_uq), ncol = b_msr$dim + 1)
    return(mat)
  }, simplify = FALSE)
  phi <- seqre_cpp(phi0, length(data),
                   b_msr$idx, b_msr$pars, b_msr$hpars,
                   s_msr$idx, s_msr$pars, s_msr$hpars,
                   eps, inc, max_it)
  out$phi <- phi
  class(out) <- c('seqre_result')
  if (keep_marg) {
    marg <- list()
    marg$args <- list(data = data, k = k, s_msr = s_msr, b_msr = b_msr,
                      keep_marg = keep_marg, burn = burn, thin = thin)
    marg$phi <- phi0
    class(marg) <- c('seqre_result')

    out <- list(full = out, marginal = marg)
  }
  return(out)
}
