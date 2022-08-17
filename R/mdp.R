#' Marginal MDP Sampler
#'
#' @param data A numeric vector of \code{n} observation values.
#' @param k The number of sampling iterations to record, will be truncated if a
#'  non-integer.
#' @param alpha The concentration parameter for the Dirichlet Process prior.
#' @param mu The mean parameter for the Normal-Inverse-Gamma prior.
#' @param tau The variance parameter for the Normal-Inverse-Gamma prior.
#' @param s The shape parameter for the Normal-Inverse-Gamma prior.
#' @param S The scale parameter for the Normal-Inverse-Gamma prior.
#' @param c The shape parameter for the Gamma prior on alpha.
#' @param C The scale parameter for the Gamma prior on alpha.
#' @param a The mean parameter for the Normal prior on mu.
#' @param A The variance parameter for the Normal prior on mu.
#' @param w The shape parameter for the Inverse-Gamma prior on tau.
#' @param W The scale parameter for the Inverse-Gamma prior on tau.
#' @param fix_a A logical value indicating whether or not to fix alpha at its
#'  initial value.
#' @param fix_m A logical value indicating whether or not to fix mu at its
#'  initial value.
#' @param fix_t A logical value indicating whether or not to fix tau at its
#'  initial value.
#' @param burn The number of initial sampling iterations to discard, will be
#'  truncated if a non-integer.
#' @param thin The number of sampling iterations to discard between records,
#'  will be truncated if a non-integer.
#'
#' @return
#' An \code{mdpolya_result} object. A list with four entries:
#'  * \code{theta}: An array of dimension \code{[k, n, 3]} encoding the
#'   component label, mean, and standard deviation for each data point for each
#'   iteration. This represents the samples from the Polya posterior
#'   distribution of the marginal MDP model.
#'  * \code{eta}: A matrix of dimension \code{[k, 5]} encoding the
#'   hyperparameter values for each iteration.
#'  * \code{args}: A list of input arguments.
#'  * \code{phi}: A list of matrices encoding the unique values from
#'   \code{theta} and associated weights for each iteration.
#' @seealso {[polya()]}
#' @references \itemize{
#'  \item Moya B., Walker S. G. (2022). Uncertainty Quantification and the
#'  Marginal MDP Model. arXiv. DOI: \doi{10.48550/arxiv.2206.08418}
#'  \item Escobar M. D., West, M. (1995) Bayesian Density Estimation and
#'  Inference Using Mixtures. Journal of the American Statistical Association.
#'  DOI: \doi{10.1080/01621459.1995.10476550}
#' }
#' @examples
#' res_mdp <- mdp(rnorm(50), 10)
#' @export
mdp <- function(data, k, alpha = 1, mu = 21, tau = 25, s = 4, S = 2,
                c = 2, C = 4, a = 21, A = 21, w = 1, W = 100,
                fix_a = FALSE, fix_m = FALSE, fix_t = FALSE,
                burn = 1000, thin = 150) {
  out <- mdp_cpp(data, k, alpha, mu, tau, s, S, c, C, a, A, w, W,
                 fix_a, fix_m, fix_t, burn, thin)
  out[[3]] <- list(data = data, k = k, alpha = alpha, mu = mu, tau = tau,
                   s = S, S = S, c = c, C = C, a = a, A = A, w = w, W = W,
                   fix_a = fix_a, fix_m = fix_m, fix_t = fix_t,
                   burn = burn, thin = thin)
  colnames(out[[1]]) <- names(data)
  dimnames(out[[1]])[[3]] <- c('z', 'mean', 'var')
  colnames(out[[2]]) <- c('alpha', 'mu', 'tau', 's', 'S')
  names(out) <- c('theta', 'eta', 'args')
  phi <- apply(out$theta, 1, function(t) {
    z <- t[, 'z']
    z_uq <- sort(unique(z))
    z_w <- table(z) / length(z)
    mat <- matrix(c(z_w, t[z_uq, 'mean'], t[z_uq, 'var']), ncol = 3)
    colnames(mat) <- c('w', 'mean', 'var')
    return(mat)
  })
  out$phi <- phi
  class(out) <- c('mdpolya_result', 'mdp')
  return(out)
}
