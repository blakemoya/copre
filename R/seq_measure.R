#' Sequence Measure for Species Sampling Models
#'
#' @param pars A list of hyperparameters used in \code{Pn} and \code{Po}.
#' @param Pn A function on a sequence length \code{n} and a number of unique
#'  values \code{k} that returns the probability of the next member in the
#'  sequence having a new value.
#' @param Po A function on a sequence length \code{n}, a number of unique values
#'  \code{k}, and the number of values equal to j \code{kj} that returns the
#'  probability of the next member in the sequence having the value j.
#' @param rnext An optional optimized generating function. If \code{NULL}, one
#'  will be generate from \code{Pn} and \code{Po}.
#'
#' @return A \code{seq_measure} object for use in the exchangeable sequence
#'  resampling scheme for mixtures.
#' @seealso {[exseqre()]}
#' @export
seq_measure <- function(pars, Pn, Po, rnext = NULL) {
  obj <- list(pars = pars, Pn = Pn, Po = Po)
  pnext <- function(z) {
    if (is.null(z) | !is.numeric(z)) {
      stop('`z` must be numeric.')
    } else if (length(z) == 0) {
      probs <- 1
      names(probs) <- 1
      return(probs)
    }
    n <- length(z)
    tab <- table(z)
    vals <- as.numeric(names(tab))
    k <- length(vals)
    probs <- c(sapply(tab, function(kj) Po(n, k, kj)), Pn(n, k))
    names(probs) <- c(vals, min(setdiff(1:(n + 1), vals)))
    return(probs)
  }
  if (is.null(rnext)) {
    rnext <- function(n, z = numeric(0)) {
      n0 <- length(z)
      zz <- numeric(n)
      for (nn in 1:n) {
        probs <- pnext(c(z, zz[0:(nn - 1)]))
        zz[nn] <- sample(as.numeric(names(probs)), 1, prob = probs)
      }
      return(zz)
    }
  }

  class(obj) <- 'seq_measure'
  obj$pnext <- pnext
  obj$rnext <- rnext
  return(obj)
}

#' A Collapsed Gnedin Process Sequence Measure.
#'
#' @param gamma The gamma parameter for the Gnedin process with xi set to 0.
#'  Bounded to \code{[0, 1]}.
#'
#' @return A \code{seq_measure} object for use in the exchangeable sequence
#'  resampling scheme for mixtures.
#' @seealso {[seq_measure(), exseqre()]}
#' @export
Sq_gnedin0 <- function(gamma) {
  Pn <- function(n, k) {
    k * (k - gamma) / (n * (gamma + n))
  }
  Po <- function(n, k, kj) {
    (gamma + n - k) / (n * (gamma + n)) * (kj + 1)
  }
  rnext <- function(n, z = numeric(0)) {
    n0 <- length(z)
    fresh <- FALSE
    if (n0 == 0) {
      z <- 1
      n0 <- 1
      fresh <- TRUE
    }
    z_uq <- unique(z)
    n_uq <- length(z_uq)
    vals <- c(sort(z_uq), setdiff(1:(n0 + n), z_uq))
    cnts <- table(factor(z, levels = vals))
    zz <- numeric(n)
    for (nn in 1:n) {
      probs <- sapply(cnts[1:n_uq], function(kj)
        Po(n0 + nn - 1, n_uq, kj)
      )
      probs <- c(probs, Pn(n0 + nn - 1, n_uq))

      zz[nn] <- sample(vals[1:(n_uq + 1)], 1, prob = probs)
      if (zz[nn] == vals[n_uq + 1]) {
        n_uq <- n_uq + 1
      }
      cnts[as.character(zz[nn])] <- cnts[as.character(zz[nn])] + 1
    }
    if (fresh) {
      return(c(z, zz[1:(n - 1)]))
    }
    return(zz)
  }
  return(seq_measure(gamma, Pn, Po, rnext))
}
