#' Base Measure for Mixture Models
#'
#' @description A structure for wrapping base measures as in Escobar and West
#' (1995).
#'
#' @param pars A list of hyperparameters used in \code{G}.
#' @param G A function which generates mixture component parameters from the
#'  conditional posterior given data \code{y}, or generates from the prior if
#'  \code{y} is empty or \code{NULL}.
#' @param q0 A marginal likelihood function on \code{y} using the values in
#'  \code{pars} (hard-coded, not as an argument).
#' @param q A likelihood function on \code{y} given \code{theta}, a mixture
#'  component parameter drawn from \code{G}.
#' @param eval An evaluation function taking \code{phi}, a list of mixture,
#'  matrices as in \code{mdp()}, \code{grd}, a grid vector, \code{f}, a
#'  character string indicating whether to calculate the gradient, density, or
#'  distribution function, and \code{nthreads}, a number of threads to utilize
#'  for parallel execution.
#'
#' @return A \code{base_measure} object for use in the exchangeable sequence
#'  resampling scheme for mixtures.
#' @seealso {[exseqre()]}
#' @references \itemize{
#'  \item Escobar M. D., West, M. (1995) Bayesian Density Estimation and
#'  Inference Using Mixtures. Journal of the American Statistical Association.
#'  DOI: \doi{10.1080/01621459.1995.10476550}
#'  }
#'  @export
base_measure <- function(pars, G, q0, q, eval) {
  obj <- list(pars = pars, G = G, q0 = q0, q = q, eval = eval)
  class(obj) <- 'base_measure'
  return(obj)
}

#' A Normal-Inverse-Gamma Base Measure for Normal Mixture Models.
#'
#' @param mu The mean parameter.
#' @param tau The variance scaling parameter.
#' @param s The primary shape parameter for the Inverse-Gamma component.
#' @param S The secondary shape parameter for the Inverse-Gamma component.
#'
#' @return A \code{base_measure} object for use in the exchangeable sequence
#'  resampling scheme for mixtures.
#' @seealso {[base_measure(), seqre()]}
#' @export
G_normgam <- function(mu, tau, s, S) {
  pars <- list(mu = mu, tau = tau, s = s, S = S)
  G <- function(y) {
    theta <- numeric(2)
    if (is.null(y) | length(y) == 0) {
      theta[1] <- rnorm(1, mu, sqrt(tau))
      theta[2] <- rgamma(1, s / 2, 2 / S)
    } else if (length(y) == 1) {
      theta[2] <- rgamma(1, (1 + s) / 2, 2 / (S + (y - mu) ^ 2 / (1 + tau)))
      x <- (mu + tau * y) / (1 + tau)
      X <- tau / (1 + tau)
      theta[1] <- rnorm(1, x, sqrt(X * theta[2]));
    } else {
      stop()
    }
    return(theta)
  }
  q0 <- function(y) {
    stopifnot(length(y) == 1)
    gamma((1 + s) / 2) / (gamma(s / 2) * sqrt(s)) *
      (1 + (y - mu) ^ 2 /
         ((1 + tau) * S)) ^ -((1 + s) / 2) /
      sqrt((1 + tau) * S / s)
  }
  q <- function(y, theta) {
    stopifnot(length(y) == 1)
    exp(-(y - theta[1]) ^ 2 / (2 * theta[2])) / sqrt(2 * theta[2])
  }
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

  return(base_measure(pars, G, q0, q, eval))
}
