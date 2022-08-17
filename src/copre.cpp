#include <RcppArmadillo.h>
#include <boost/math/special_functions/erf.hpp>
#if defined(_OPENMP)
  #include <omp.h>
  // [[Rcpp::plugins(openmp)]]
#endif

static const double NC_NORM = 0.3989422804014327;

double normcdf(const double x) {
  double value = 0.5 * erfc(-x * M_SQRT1_2);
  return value;
}

double normcdfinv(const double x) {
  double value = -M_SQRT2 * boost::math::erfc_inv(2 * x);
  return value;
}

double normpdf(const double x) {
  double value = NC_NORM * exp(-0.5 * x * x);
  return value;
}

double bivnormpdf(const double x, const double y, const double rho) {
  const double NC_BIVNORM = NC_NORM * NC_NORM / sqrt(1 - rho * rho);
  const double z = x * x + y * y - 2.0 * rho * x * y;
  double value =  NC_BIVNORM * exp(-0.5 * z / (1 - rho * rho));
  return value;
}

double dbgc(const double u, const double v, const double rho) {
  const double x = normcdfinv(u);
  const double y = normcdfinv(v);
  double value = bivnormpdf(x, y, rho) / (normpdf(x) * normpdf(y));
  return value;
}

double H(const double u, const double v, const double rho) {
  const double x = normcdfinv(u);
  const double y = normcdfinv(v);
  double value = normcdf((x - rho * y) / sqrt(1 - rho * rho));
  return value;
}

//[[Rcpp::export]]
arma::mat copre_cpp(arma::mat y_perm, arma::vec alpha, double rho,
                    arma::uword N, arma::uword k, arma::mat &P, arma::vec grd,
                    arma::uword nthreads) {
  arma::uword n = y_perm.n_rows;
  arma::mat u(k, N, arma::fill::randu);
  #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthreads)
  #endif
  for (arma::uword kk = 0; kk < k; kk++) {
    arma::vec y = y_perm.col(kk);
    arma::vec Pnn(n);
    double obs;
    for (arma::uword nn = 0; nn < n; nn++) {
      arma::interp1(grd, P.col(kk), y, Pnn);
      obs = Pnn(nn);
      for (arma::uword m = 0; m < grd.n_elem; m++) {
        double past = (1 - alpha(nn)) * P(m, kk);
        double pres = alpha(nn) * H(P(m, kk), obs, rho);
        P(m, kk) = past + pres;
      }
    }
  }
  #if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthreads) collapse(2)
  #endif
  for (arma::uword kk = 0; kk < k; kk++) {
    for (arma::uword m = 0; m < grd.n_elem; m++) {
      double obs;
      for (arma::uword NN = 0; NN < N; NN++) {
        obs = u(kk, NN);
        double past = (1 - alpha(n + NN)) * P(m, kk);
        double pres = alpha(n + NN) * H(P(m, kk), obs, rho);
        P(m, kk) = past + pres;
      }
    }
  }
  return P.t();
}
