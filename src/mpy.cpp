#include "base.h"

#if defined(_OPENMP)
  #include <omp.h>
#endif

// [[Rcpp::export]]
Rcpp::List mpy_cpp(arma::vec y, arma::uword k,
                   double d, double alpha, double mu, double tau,
                   double s, double S, double c, double C, double a, double A,
                   double w, double W, bool fix_a, bool fix_m, bool fix_t,
                   arma::uword burn, arma::uword thin) {
  arma::uword n = y.n_elem;
  arma::cube recs(k, n, 3);
  arma::mat hprecs(k, 5);
  int rec_cnt = 0;
  arma::uvec z(n, arma::fill::ones);
  z = cumsum(z) - 1;
  arma::uvec seq = z;
  arma::mat theta(n, 2);
  for (arma::uword nn = 0; nn < n; nn++) {
    theta.row(nn) = rG(y(nn), s, S, mu, tau).t();
  }
  for (arma::uword iter = 0; iter < burn + thin * k; iter++) {
    arma::uvec z_uq = unique(z);
    arma::mat theta_uq = theta.rows(z_uq);
    arma::uword n_uq = z_uq.n_elem;
    if (!fix_m) {
      double x = A / (A + tau * sum(1 / theta_uq.col(1)));
      double a_post = (1 - x) * a + x * sum(theta_uq.col(1))
        * sum(theta_uq.col(0) / theta_uq.col(1));
      double A_post = x * tau * sum(theta_uq.col(1));
      mu = R::rnorm(a_post, sqrt(A_post));
    }
    if (!fix_t) {
      double w_post = w + n_uq;
      double W_post = W + sum(pow(theta_uq.col(0) - mu, 2.0) / theta_uq.col(1));
      tau = 1.0 / R::rgamma(w_post / 2.0, 2.0 / W_post);
    }
    for (arma::uword i = 0; i < n; i++) {
      double q0 = q(y(i), d, alpha, s, S, mu, tau, n - 1);
      arma::vec qj(n, arma::fill::zeros);
      for (arma::uword j = 0; j < n; j++) {
        if (j == i) {
          qj(j) = q0;
        } else {
          qj(j) = (1 - d) * q(y(i), theta(z(j), 0), theta(z(j), 1));
        }
      }
      qj = qj / sum(qj);
      arma::uword z_prop = sample(seq, 1, true, qj)(0);
      if (z_prop == i) {
        theta.row(i) = rG(y(i), s, S, mu, tau).t();
        z(i) = z_prop;
      } else {
        z(i) = z(z_prop);
      }
    }
    if (!fix_a) {
      double eta = R::rbeta(alpha + 1, n);
      double u = R::runif(0, 1);
      if (u < (c + n_uq - 1) / (n * (C - log(eta)) + c + n_uq - 1)) {
        alpha = R::rgamma(c + n_uq, 1.0 / (C - log(eta)));
      } else {
        alpha = R::rgamma(c + n_uq - 1, 1.0 / (C - log(eta)));
      }
    }
    if ((iter >= burn) && (iter % thin == 0)) {
      recs.slice(0).row(rec_cnt) = arma::conv_to<arma::vec>::from(z + 1).t();
      recs.slice(1).row(rec_cnt) = theta.col(0).t();
      recs.slice(2).row(rec_cnt) = theta.col(1).t();
      hprecs.row(rec_cnt) = arma::vec({alpha, mu, tau, s, S}).t();
      rec_cnt++;
    }
  }
  return Rcpp::List::create(recs, hprecs);
}
