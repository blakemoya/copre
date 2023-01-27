#' Base Measure for Mixture Models
#'
#' @description A structure for wrapping base measures as in Escobar and West
#' (1995).
#'
#' @param idx A unique index for the base measure.
#' @param dim A dimension for the support of the base measure.
#' @param pars A list of parameters used to generate mixture components.
#' @param hpars A list of hyperparameters used to generate `pars`.
#' @param eval An evaluation function taking `phi`, a list of mixture parameter
#'   matrices, `grd`, a grid vector, `f`, a character string indicating whether
#'   to calculate the gradient, density, or distribution function, and
#'   `nthreads`, a number of threads to utilize for parallel execution.
#'
#' @return A `base_measure` object for use in the sequence resampling scheme for
#'   mixtures.
#' @seealso {[seqre()]}
#' @references \itemize{
#'  \item Escobar M. D., West, M. (1995) Bayesian Density Estimation and
#'  Inference Using Mixtures. Journal of the American Statistical Association.
#'  DOI: \doi{10.1080/01621459.1995.10476550}
#'  }
base_measure <- function(idx, dim, pars, hpars, eval) {
  obj <- list(idx = idx, dim = dim, pars = pars, hpars = hpars, eval = eval)
  class(obj) <- 'base_measure'
  return(obj)
}

#' Normal-Inverse-Gamma Base Measure for Location-Scale Normal Mixture Models.
#'
#' @param mu The mean parameter.
#' @param tau The variance scaling parameter.
#' @param s The primary shape parameter for the Inverse-Gamma component.
#' @param S The secondary shape parameter for the Inverse-Gamma component.
#' @param a The prior mean parameter for `mu`.
#' @param A The prior variance for `mu`.
#' @param fix_m A logical value indicating whether or not `mu` should be fixed.
#' @param w The prior primary shape parameter for `tau`.
#' @param W The prior secondary shape parameter for `tau`.
#' @param fix_t A logical value indicating whether or not `tau` should be fixed.
#'
#' @return A `base_measure` object for use in the sequence resampling scheme for
#'   mixtures.
#' @seealso {[base_measure(), seqre()]}
#' @export
G_normls <- function(mu = 0, tau = 1, s = 1, S = 1,
                     a = NULL, A = NULL, w = NULL, W = NULL) {
  if (!(is.null(a) | is.null(A))) {
    fix_m <- FALSE
    mu <- rnorm(1, a, sqrt(A))
  } else {
    fix_m <- TRUE
  }
  if (!(is.null(w) | is.null(W))) {
    fix_t <- FALSE
    tau <- rgamma(1, w, W)
  } else {
    fix_t <- TRUE
  }
  pars <- c(mu = mu, tau = tau, s = s, S = S)
  hpars <- c(a = a, A = A, fix_m = fix_m, w = w, W = W, fix_t = fix_t)
  eval <- function(phi, grd, f = 'density', nthreads = 1) {
    out <- sapply(phi, function(p) {
      if (f == 'gradient') {
        stop('`gradient` not implemented for `G_normgam`.')
      } else if (f == 'density') {
        func <- function(x) {
          sum(p[, 1] * dnorm(x, p[, 2], sqrt(p[, 3])))
        }
      } else if (f == 'distribution') {
        func <- function(x) {
          sum(p[, 1] * pnorm(x, p[, 2], sqrt(p[, 3])))
        }
      }
      return(sapply(grd, func))
    })
    return(t(out))
  }

  return(base_measure(0, 2, pars, hpars, eval))
}
