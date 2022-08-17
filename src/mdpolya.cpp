#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
#include <boost/math/distributions.hpp>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
#if defined(_OPENMP)
#include <omp.h>
// // [[Rcpp::plugins(openmp)]] // This causes problems on MacOS
#endif

using namespace arma;


arma::vec rG(double s, double S, double mu, double tau) {
  arma::vec theta(2);
  theta(1) = 1.0 / R::rgamma(s / 2.0, 2.0 / S);
  theta(0) = R::rnorm(mu, sqrt(tau));
  return theta;
}

arma::vec rG(double y, double s, double S, double mu, double tau) {
  arma::vec theta(2);
  theta(1) = 1.0 / R::rgamma((1.0 + s) / 2.0,
        2.0 / (S + pow(y - mu, 2.0) / (1.0 + tau)));
  double x = (mu + tau * y) / (1 + tau);
  double X = tau / (1 + tau);
  theta(0) = R::rnorm(x, sqrt(X * theta(1)));
  return theta;
}

double q(double y, double alpha, double s, double S, double mu, double tau) {
  double q0 = alpha * tgamma((1.0 + s) / 2.0) / (tgamma(s / 2.0) * sqrt(s));
  q0 *= pow(1.0 + pow(y - mu, 2.0) / ((1 + tau) * S), -(1 + s) / 2.0);
  q0 /= sqrt((1 + tau) * S / s);
  return q0;
}

double q(double y, double mu, double v) {
  return exp(-pow(y - mu, 2.0) / (2 * v)) / sqrt(2 * v);
}

// [[Rcpp::export]]
Rcpp::List mdp_cpp(arma::vec y, arma::uword k,
                   double alpha, double mu, double tau, double s, double S,
                   double c, double C, double a, double A, double w, double W,
                   bool fix_a, bool fix_m, bool fix_t,
                   arma::uword burn, arma::uword thin) {
  arma::uword n = y.n_elem;
  arma::cube recs(k, n, 3);
  arma::mat hprecs(k, 5);
  int rec_cnt = 0;
  arma::uvec z(n, fill::ones);
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
      double q0 = q(y(i), alpha, s, S, mu, tau);
      arma::vec qj(n, fill::zeros);
      for (arma::uword j = 0; j < n; j++) {
        if (j == i) {
          qj(j) = q0;
        } else {
          qj(j) = q(y(i), theta(z(j), 0), theta(z(j), 1));
        }
      }
      qj = qj / sum(qj);
      arma::uword z_prop = Rcpp::RcppArmadillo::sample(seq, 1, true, qj)(0);
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
      recs.slice(0).row(rec_cnt) = conv_to<arma::vec>::from(z + 1).t();
      recs.slice(1).row(rec_cnt) = theta.col(0).t();
      recs.slice(2).row(rec_cnt) = theta.col(1).t();
      hprecs.row(rec_cnt) = arma::vec({alpha, mu, tau, s, S}).t();
      rec_cnt++;
    }
  }
  return Rcpp::List::create(recs, hprecs);
}

// [[Rcpp::export]]
Rcpp::List polyaurn_cpp(arma::cube theta, arma::mat hyper, double eps, double ups,
                        arma::uword nthreads) {
  arma::uword k = hyper.n_rows;
  arma::uword n = theta.n_cols;
  std::vector<arma::mat> out(k);
  #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthreads)
  #endif
  for (arma::uword kk = 0; kk < k; kk++) {
    double alpha = hyper(kk, 0);
    arma::uword m = 1.0 + R::qpois(1 - ups, -(alpha + n) * std::log(eps),
                             true, false);
    arma::mat phi(m, 3, fill::zeros);
    arma::vec v(m);
    for (arma::uword mm = 0; mm < m; mm++) {
      v(mm) = R::rbeta(1, alpha + n);
    }
    phi.col(0) = v;
    v = cumprod(1 - v);
    arma::rowvec phi_ex(3, fill::zeros);
    phi = join_cols(phi, phi_ex);
    arma::uvec z_track(n, fill::value(m + 1));
    for (arma::uword mm = 0; mm <= m; mm++) {
      if (mm > 0) {
        if (mm == m) {
          phi(mm, 0) = 1.0 - sum(phi.col(0));
        } else {
          phi(mm, 0) *= v(mm - 1);
        }
      }
      double u = R::runif(0, 1) * (alpha + n);
      if (u <= alpha) {
        phi(mm, 1) = R::rnorm(hyper(kk, 1), std::sqrt(hyper(kk, 2)));
        phi(mm, 2) = R::rgamma(hyper(kk, 3) / 2.0, 2.0 / hyper(kk, 4));
      } else {
        int i = std::floor(R::runif(0, 1) * n);
        int z = theta(kk, i, 0) - 1;
        if (z_track(z) > m) {
          phi(mm, 1) = theta(kk, z, 1);
          phi(mm, 2) = theta(kk, z, 2);
          z_track(z) = mm;
        } else {
          arma::uword mmm = z_track(z);
          phi(mmm, 0) += phi(mm, 0);
          phi(mm, 0) = 0.0;
          phi(mm, 1) = datum::nan;
          phi(mm, 2) = datum::nan;
        }
      }
    }
    out[kk] = phi;
  }
  Rcpp::List out_;
  for (arma::uword kk = 0; kk < k; kk++) {
    out_.push_back(out[kk]);
  }
  return out_;
}

// [[Rcpp::export]]
Rcpp::List polyaurngen_cpp(arma::umat theta, arma::vec alpha, double eps, double ups,
                           arma::uword nthreads) {
  arma::uword k = theta.n_rows;
  arma::uword n = theta.n_cols;
  std::vector<arma::mat> out(k);
  #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthreads)
  #endif
  for (arma::uword kk = 0; kk < k; kk++) {
    arma::uword m = 1.0 + R::qpois(1 - ups, -(alpha(kk) + n) * std::log(eps),
                             true, false);
    arma::vec v(m);
    arma::mat phi(m, 2, fill::zeros);
    for (arma::uword mm = 0; mm < m; mm++) {
      v(mm) = R::rbeta(1, alpha(kk) + n);
    }
    phi.col(0) = v;
    v = cumprod(1 - v);
    arma::rowvec phi_ex(2, fill::zeros);
    phi = join_cols(phi, phi_ex);
    arma::uword z_max = max(theta.row(kk));
    arma::uvec z_track(z_max, fill::value(m + 1));
    for (arma::uword mm = 0; mm <= m; mm++) {
      if (mm > 0) {
        if (mm == m) {
          phi(mm, 0) = 1.0 - sum(phi.col(0));
        } else {
          phi(mm, 0) *= v(mm - 1);
        }
      }
      double u = R::runif(0, 1) * (alpha(kk) + n);
      if (u <= alpha(kk)) {
        phi(mm, 1) = z_max;
        z_max++;
      } else {
        int i = std::floor(R::runif(0, 1) * n);
        int z = theta(kk, i) - 1;
        if (z_track(z) > m) {
          phi(mm, 1) = z;
          z_track(z) = mm;
        } else {
          arma::uword mmm = z_track(z);
          phi(mmm, 0) += phi(mm, 0);
          phi(mm, 0) = 0.0;
          phi(mm, 1) = datum::nan;
        }
      }
    }
    out[kk] = phi;
  }
  Rcpp::List out_;
  for (arma::uword kk = 0; kk < k; kk++) {
    out_.push_back(out[kk]);
  }
  return out_;
}

double grad(double x, double mean, double var) {
  return R::dnorm(x, mean, var, false) * (mean - x) / var;
}

// [[Rcpp::export]]
arma::mat evalmdpolya_cpp(Rcpp::List &mdpolya, arma::vec &x, arma::uword f, arma::uword nthreads) {
  arma::mat out(mdpolya.length(), x.n_elem, fill::zeros);
  #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthreads)
  #endif
  for (arma::uword kk = 0; kk < mdpolya.length(); kk++) {
    arma::mat thetak = mdpolya[kk];
    for (arma::uword xx = 0; xx < x.n_elem; xx++) {
      for (arma::uword zz = 0; zz < thetak.n_rows; zz++) {
        if (f == 0) {
          out(kk, xx) += thetak(zz, 0) * R::pnorm(x(xx), thetak(zz, 1),
              std::sqrt(thetak(zz, 2)), true, false);
        } else if (f == 1) {
          out(kk, xx) += thetak(zz, 0) * R::dnorm(x(xx), thetak(zz, 1),
              std::sqrt(thetak(zz, 2)), false);
        } else if (f == 2) {
          out(kk, xx) += thetak(zz, 0) * grad(x(xx), thetak(zz, 1),
              thetak(zz, 2));
        }
      }
    }
  }
  return(out);
}
