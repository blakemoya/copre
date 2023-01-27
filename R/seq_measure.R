#' Sequence Measure for Species Sampling Models
#'
#' @param idx A unique index for the sequence measure.
#' @param pars A list of parameters used in `Pn` and `Po` to generate a
#'   sequence.
#' @param hpars A list of hyperparameters used to generate `pars`.
#' @param Pn A function on a sequence length `n` and a number of unique values
#'   `k` that returns the probability of the next member in the sequence having
#'   a new value.
#' @param Po A function on a sequence length `n`, a number of unique values `k`,
#'   and the number of values equal to `j`, `kj`, that returns the probability
#'   of the next member in the sequence having the value `j`.
#'
#' @return A `seq_measure` object for use in the exchangeable sequence
#'   resampling scheme for mixtures.
#' @seealso [seqre()]
seq_measure <- function(idx, pars, hpars, Pn, Po) {
  obj <- list(idx = idx, pars = pars, hpars = hpars, Pn = Pn, Po = Po)
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

  class(obj) <- 'seq_measure'
  obj$pnext <- pnext
  obj$rnext <- rnext
  return(obj)
}

#' Dirichlet Sequence Measure.
#'
#' @param alpha The concentration parameter for the Dirichlet process. Must be
#'   greater than 0.
#' @param c The prior primary shape parameter for `alpha`.
#' @param C The prior secondary shape parameter for `alpha`.
#' @param fix_m A logical value indicating whether or not `alpha` should be
#'   fixed.
#'
#' @return A `seq_measure` object for use in the exchangeable sequence
#'   resampling scheme for mixtures.
#' @seealso [seq_measure()], [seqre()]
#' @export
Sq_dirichlet <- function(alpha = 1, c = NULL, C = NULL) {
  if (!(is.null(c) | is.null(C))) {
    fix_a <- FALSE
    alpha <- rgamma(c, C)
  } else {
    fix_a <- TRUE
  }
  if (alpha <= 0) {
    stop('`alpha` must be positive.')
  }
  Pn <- function(n, k) {
    alpha / (alpha + n)
  }
  Po <- function(n, k, kj) {
    kj / (alpha + n)
  }
  return(seq_measure(0, alpha, c(c = c, C = C, fix_a = fix_a), Pn, Po))
}

#' Pitman-Yor Sequence Measure.
#'
#' @param d The discount parameter for the Pitman-Yor process. Must be less than
#'   1.
#' @param alpha The concentration parameter for the Pitman-Yor process. Must be
#'   greater than -`sigma` if `sigma` is in [0, 1), else ignored.
#' @param m A positive integer used to set `theta = m * abs(sigma)` if `sigma`
#'   is negative.
#'
#' @return A \code{seq_measure} object for use in the exchangeable sequence
#'  resampling scheme for mixtures.
#' @seealso [seq_measure()], [seqre()]
#' @export
Sq_pitmanyor <- function(d, alpha = 1, m = 1L) {
  if (d < 0) {
    if (m < 0 | !is.integer(m)) {
      stop('`m` must be a positive integer. Remember to add `L` after the
           numeric value.')
    }
    alpha <- m * abs(d)
  }
  Pn <- function(n, k) {
    (alpha + d * k) / (alpha + n)
  }
  Po <- function(n, k, kj) {
    (kj - d)  / (alpha + n)
  }
  return(seq_measure(1, c(d = d, alpha = alpha), 0, Pn, Po))
}

#' Collapsed Gnedin Process Sequence Measure.
#'
#' @param gamma The gamma parameter for the Gnedin process with xi set to 0.
#'   Bounded to `[0, 1]`.
#'
#' @return A `seq_measure` object for use in the exchangeable sequence
#'   resampling scheme for mixtures.
#' @seealso [seq_measure()], [seqre()]
#' @export
Sq_gnedin0 <- function(gamma) {
  Pn <- function(n, k) {
    k * (k - gamma) / (n * (gamma + n))
  }
  Po <- function(n, k, kj) {
    (gamma + n - k) / (n * (gamma + n)) * (kj + 1)
  }
  return(seq_measure(2, c(gamma = gamma), 0, Pn, Po))
}
