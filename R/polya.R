#' Polya Completion for Marginal MDP Samples
#'
#' @param res_mdp Samples from a marginal MDP model.
#' @param epsilon The desired maximum weight associated with the final remainder
#'   component.
#' @param upsilon The portion of samples which do not meet the desired epsilon.
#' @param nthreads UNSTABLE: The number of parallel threads to launch with
#'   OpenMP, not recommended due to induced instability.
#'
#' @return If \code{res_mdp} was an \code{mdpolya_result} object, returns
#'   another \code{mdpolya_result}object with \code{phi}, \code{eta} and
#'   \code{args} entries as in [mdp()]. If \code{res_mdp} was a
#'   \code{dirichletprocess} object, returns another \code{dirichletprocess}
#'   object with new components and altered weights.
#'
#' @references Moya B., Walker S. G. (2022). Uncertainty Quantification and the
#'  Marginal MDP Model. arXiv. DOI: \doi{10.48550/arxiv.2206.08418}
#' @examples
#' res_mdp <- mdp(rnorm(50), 10)
#' res_pol <- polya(res_mdp, nthreads = 1)
#' @export
polya <- function(res_mdp, epsilon = 0.01, upsilon = 0.01,
                  nthreads = 1) {
  UseMethod('polya')
}

#' @describeIn polya Polya extension to a \code{mdpolya_result} object.
#' @export
polya.mdp <- function(res_mdp, epsilon = 0.01, upsilon = 0.01,
                      nthreads = 1) {
  stopifnot(any(class(res_mdp) == 'mdp') &
              epsilon > 0 & epsilon < 1 &
              upsilon > 0 & upsilon < 1)
  out <- list(phi = NULL, eta = res_mdp$eta, args = res_mdp$args);
  out$phi <- polyaurn_cpp(res_mdp$theta, res_mdp$eta,
                          epsilon, upsilon, nthreads)
  sapply(1:length(out$phi), function(i) {
    out$phi[[i]] <<- out$phi[[i]][is.finite(out$phi[[i]][, 2]), , drop = FALSE];
    if (!is.matrix(out$phi[[i]])) { # This may be unnecessary now due to drop
      out$phi[[i]] <<- matrix(out$phi[[i]], nrow = 1)
    }
    colnames(out$phi[[i]]) <<- c('w', 'mean', 'var')
  })
  class(out) <- c('mdpolya_result', 'polya')
  return(out)
}

#' @describeIn polya Polya extension to a \code{dirichletprocess} object.
#' @export
polya.dirichletprocess <- function(res_mdp, epsilon = 0.01, upsilon = 0.01,
                                   nthreads = 1) {
  phi <- polyaurngen_cpp(t(sapply(res_mdp$labelsChain, identity)),
                         res_mdp$alphaChain, epsilon, upsilon, nthreads)
  sapply(1:length(phi), function(i) {
    phi[[i]] <<- phi[[i]][is.finite(phi[[i]][, 2]), ];
    if (!is.matrix(phi[[i]])) {
      phi[[i]] <<- matrix(phi[[i]], nrow = 1)
    }
    colnames(phi[[i]]) <<- c('w', 'z')
    phi[[i]][, 'z'] <<- phi[[i]][, 'z'] + 1
  })
  res_mdp$weightsChain <- lapply(phi, function(p) p[, 'w'])
  res_mdp$clusterParametersChain <- lapply(1:length(phi), function(i) {
    ns <- phi[[i]][, 'z'] > max(res_mdp$labelsChain[[i]])
    if(sum(ns) != 0) {
      pars <- dirichletprocess::PriorDraw(res_mdp$mixingDistribution, sum(ns))
    } else {
      pars <- list(array(dim = c(0, 0, 0)))
    }
    out <- lapply(1:length(res_mdp$clusterParametersChain[[i]]), function(j) {
      if (!all(!ns)) {
        if (prod(dim(pars[[j]])) != 0) {
          abind::abind(res_mdp$clusterParametersChain[[i]][[j]][, , phi[[i]][, 'z'][!ns],
                                                                drop = FALSE],
                       pars[[j]])
        } else {
          res_mdp$clusterParametersChain[[i]][[j]][, , which(!ns),drop = FALSE]
        }
      } else {
        res_mdp$clusterParametersChain[[i]][[j]][, , which(!ns),drop = FALSE]
      }
    })
    return(out)
  })
  return(res_mdp)
}
