#include <RcppArmadillo.h>
#include <math.h>
#include <boost/math/distributions.hpp>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
#if defined(_OPENMP)
  #include <omp.h>
  // // [[Rcpp::plugins(openmp)]] // This causes problems on MacOS
#endif


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
    arma::mat phi(m, 2, arma::fill::zeros);
    for (arma::uword mm = 0; mm < m; mm++) {
      v(mm) = R::rbeta(1, alpha(kk) + n);
    }
    phi.col(0) = v;
    v = cumprod(1 - v);
    arma::rowvec phi_ex(2, arma::fill::zeros);
    phi = join_cols(phi, phi_ex);
    arma::uword z_max = max(theta.row(kk));
    arma::uvec z_track(z_max, arma::fill::value(m + 1));
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
          phi(mm, 1) = arma::datum::nan;
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
  arma::mat out(mdpolya.length(), x.n_elem, arma::fill::zeros);
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
