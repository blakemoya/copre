#' Sequence Resampling
#'
#' @description A function that samples predictive distributions for univariate
#'  continuous data using exchangeable predictive extension.
#'
#' @param data The data from which to sample predictive distributions.
#' @param k The number of predictive samples to draw.
#' @param s_msr A \code{seq_measure} object.
#' @param b_msr A univariate \code{base_measure} object.
#' @param burn The number of initial sampling iterations to discard, will be
#'  truncated if a non-integer.
#' @param thin The number of sampling iterations to discard between records,
#'  will be truncated if a non-integer.
#'
#' @return A \code{unimplemented} object.
#' @seealso [seq_measure()], [base_measure()]
#' @export
seqre <- function(data, k, s_msr, b_msr, marg = FALSE, burn = 1000,
                  thin = 150) {
  n <- length(data)
  recs <- array(NA, dim = c(k, n, 3))
  stopifnot('seq_measure' %in% class(s_msr))
  stopifnot('base_measure' %in% class(b_msr))

  z <- s_msr$rnext(n)
  theta <- t(sapply(1:n, function(nn) b_msr$G(NULL)))
  th_dim <- ncol(theta)
  kk <- 1

  for (iter in 1:(burn + thin * k)) {
    for (i in 1:n) {
      zi <- z[1:n != i]
      z_uq <- sort(unique(zi))
      k1 <- length(z_uq)
      q0 <- b_msr$q0(data[i]) * s_msr$Pn(n - 1, k1)
      qj <- sapply(z_uq, function(zz) {
        b_msr$q(data[i], theta[zz, ]) * s_msr$Po(n - 1, k1, sum(zi == zz))
      })
      qj <- c(qj, q0)
      z_new <- min(setdiff(1:n, z_uq))
      z_prop <- sample(c(z_uq, z_new), 1, prob = qj)
      if (z_prop == z_new) {
        theta[z_new, ] <- b_msr$G(data[i])
        z[i] <- z_new
      } else {
        z[i] <- z_prop
      }
    }
    if ((iter > burn) & ((iter - burn) %% thin == 0)) {
      recs[kk, , 1] <- z
      recs[kk, , 2:dim(recs)[[3]]] <- theta
      kk <- kk + 1
    }
  }

  out <- list()
  out$theta <- recs
  out$args <- list(data = data, k = k, s_msr = s_msr, b_msr = b_msr,
                   marg = marg, burn = burn, thin = thin)
  phi <- apply(recs, 1, function(t) {
    z <- t[, 1]
    z_uq <- sort(unique(z))
    z_w <- table(z) / length(z)
    mat <- matrix(c(z_w, t[z_uq, 2:(th_dim + 1)]),
                  nrow = length(z_uq), ncol = th_dim + 1)
    return(mat)
  }, simplify = FALSE)
  if(!marg) {
    N <- 1000
    epsilon <- 0.001
    phi <- lapply(phi, function(p) {
      z <- rep(1:nrow(p), p[, 1] * n)
      z_uq <- sort(unique(z))
      z_w <- table(z) / length(z)
      w_max <- max(z_w)
      w_diff <- 1
      while(w_diff > epsilon) {
        z <- c(z, s_msr$rnext(N, z))
        z_w <- table(z) / length(z)
        w_max_ <- max(z_w)
        w_diff <- abs(w_max - w_max_)
        w_max <- w_max_
      }
      mat <- matrix(NA, nrow = length(z_w), ncol = th_dim + 1)
      m <- nrow(mat) - nrow(p)
      if (m == 0) {
        p[, 1] <- z_w
        return(p)
      }
      mat[, 1] <- z_w
      mat[, 2:(th_dim + 1)] <- rbind(p[, 2:(th_dim + 1)],
                                     t(sapply(1:m, function(mm) b_msr$G(NULL))))
      return(mat)
    })
  }
  out$phi <- phi
  class(out) <- c('seqre_result')
  return(out)
}
