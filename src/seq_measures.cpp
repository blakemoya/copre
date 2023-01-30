#include "seq_measures.h"

arma::uvec setdiff(arma::uvec x, const arma::uvec& y) {
  for (size_t j = 0; j < y.n_elem; j++) {
    arma::uword q1 = arma::conv_to<arma::uword>::from(arma::find(x == y[j]));
    x.shed_row(q1);
  }
  return x;
}

arma::uvec rnext(seq_measure& s, arma::uword n, arma::uvec cnts, bool tab) {
  arma::uword n_old = arma::sum(cnts);
  arma::uvec nums(n, arma::fill::ones);
  nums = cumsum(nums) - 1;
  arma::uvec z(n, arma::fill::zeros);
  arma::uword zz = arma::max(arma::find(cnts > 0)) + 1;
  arma::vec probs(n);
  for (arma::uword nn = 0; nn < n; nn++) {
    for (arma::uword kk = 0; kk < zz; kk++) {
      probs(kk) = s.p_old(nn + n_old, zz - 1, cnts(kk));
    }
    probs(zz) = s.p_new(nn + n_old, zz - 1);
    arma::uvec z_next_ = sample(nums, 1, false, probs);
    arma::uword z_next = z_next_(0);
    if (z_next == zz) {
      zz++;
      cnts = join_cols(cnts, arma::uvec({0}));
    }
    cnts(z_next)++;
    z(nn) = z_next;
  }
  if (tab) {
    arma::uword maxk = max(find(cnts > 0));
    return cnts(arma::span(0, maxk));
  } else {
    return z;
  }
}

arma::uvec rnext(seq_measure& s, arma::uword n, bool tab) {
  arma::uvec nums(n, arma::fill::ones);
  nums = cumsum(nums) - 1;
  arma::uvec z(n, arma::fill::zeros);
  arma::uvec cnts(n, arma::fill::zeros);
  z(0) = 0;
  cnts(0)++;
  arma::uword zz = 1;
  arma::vec probs(n);
  for (arma::uword nn = 1; nn < n; nn++) {
    for (arma::uword kk = 0; kk < zz; kk++) {
      probs(kk) = s.p_old(nn, zz - 1, cnts(kk));
    }
    probs(zz) = s.p_new(nn, zz - 1);
    arma::uvec z_next_ = sample(nums, 1, false, probs);
    arma::uword z_next = z_next_(0);
    if (z_next == zz) {
      zz++;
    }
    cnts(z_next)++;
    z(nn) = z_next;
  }
  if (tab) {
    arma::uword maxk = max(find(cnts > 0));
    return cnts(arma::span(0, maxk));
  } else {
    return z;
  }
}

bool sq_gnedin0::validate() {
  return (gamma <= 0) & (gamma <= 1);
}

double sq_gnedin0::p_new(arma::uword n, arma::uword k) {
  return k * (k - gamma) / (n * (gamma + n));
}

double sq_gnedin0::p_old(arma::uword n, arma::uword k, arma::uword kj) {
  return (gamma + n - k) / (n * (gamma + n)) * (kj + 1.0);
}

arma::vec sq_gnedin0::pars() {
  arma::vec out(1);
  out(0) = gamma;
  return out;
}

arma::vec sq_gnedin0::hpars() {
  arma::vec out;
  return out;
}

arma::vec sq_gnedin0::update(arma::uword n, arma::uvec cnts) {
  arma::vec out(1);
  out(0) = gamma;
  return out;
}

bool sq_pitmanyor::validate() {
  if (d < 0) {
    alpha = m * std::abs(d);
    return true;
  } else {
    return (d < 1) & (alpha > -d);
  }
}

double sq_pitmanyor::p_new(arma::uword n, arma::uword k) {
  return (alpha + d * k) / (alpha + n);
}

double sq_pitmanyor::p_old(arma::uword n, arma::uword k, arma::uword kj) {
  return (kj - d) / (alpha + n);
}

arma::vec sq_pitmanyor::pars() {
  arma::vec out(2);
  out(0) = d;
  out(1) = alpha;
  return out;
}

arma::vec sq_pitmanyor::hpars() {
  arma::vec out;
  return out;
}

arma::vec sq_pitmanyor::update(arma::uword n, arma::uvec cnts) {
  arma::vec out(2);
  out(0) = d;
  out(1) = alpha;
  return out;
}

bool sq_dirichlet::validate() {
  return (alpha > 0);
}

double sq_dirichlet::p_new(arma::uword n, arma::uword k) {
  return alpha / (alpha + n);
}

double sq_dirichlet::p_old(arma::uword n, arma::uword k, arma::uword kj) {
  return kj / (alpha + n);
}

arma::vec sq_dirichlet::pars() {
  arma::vec out(1);
  out(0) = alpha;
  return out;
}

arma::vec sq_dirichlet::hpars() {
  arma::vec out(2);
  out(0) = c;
  out(1) = C;
  return out;
}

arma::vec sq_dirichlet::update(arma::uword n, arma::uvec cnts) {
  if (!fix_a) {
    arma::uword k = cnts.n_elem;
    double eta = R::rbeta(alpha + 1, n);
    double u = R::runif(0, 1);
    if (u < (c + k - 1) / (n * (C - log(eta)) + c + k - 1)) {
      alpha = R::rgamma(c + cnts.n_elem, 1.0 / (C - log(eta)));
    } else {
      alpha = R::rgamma(c + k - 1, 1.0 / (C - log(eta)));
    }
  }
  arma::vec out(1);
  out(0) = alpha;
  return out;
}
