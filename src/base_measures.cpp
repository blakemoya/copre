#include "base_measures.h"

bool g_normls::validate() {
  return (tau > 0) & (s > 0) &  (S > 0) & (A > 0) & (w > 0) & (W > 0);
}

arma::vec g_normls::rtheta() {
  arma::vec theta(2);
  theta(1) = 1.0 / R::rgamma(s / 2.0, 2.0 / S);
  theta(0) = R::rnorm(mu, sqrt(tau));
  return theta;
}

arma::vec g_normls::rtheta(double y) {
  arma::vec theta(2);
  theta(1) = 1.0 / R::rgamma((1.0 + s) / 2.0,
        2.0 / (S + pow(y - mu, 2.0) / (1.0 + tau)));
  double x = (mu + tau * y) / (1 + tau);
  double X = tau / (1 + tau);
  theta(0) = R::rnorm(x, sqrt(X * theta(1)));
  return theta;
}

double g_normls::q(double y) {
  double q0 = tgamma((1.0 + s) / 2.0) / (tgamma(s / 2.0) * sqrt(s));
  q0 *= pow(1.0 + pow(y - mu, 2.0) / ((1 + tau) * S), -(1 + s) / 2.0);
  q0 /= sqrt((1 + tau) * S / s);
  return q0;
}

double g_normls::q(double y, arma::vec theta) {
  return exp(-pow(y - theta(0), 2.0) / (2 * theta(1))) / sqrt(2 * theta(1));
}

arma::vec g_normls::pars() {
  arma::vec out(4);
  out(0) = mu;
  out(1) = tau;
  out(2) = s;
  out(3) = S;
  return out;
}

arma::vec g_normls::hpars() {
  arma::vec out(4);
  out(0) = a;
  out(1) = A;
  out(2) = w;
  out(3) = W;
  return out;
}

arma::vec g_normls::update(arma::mat theta) {
  if (!fix_m) {
    double x = A / (A + tau * sum(1 / theta.col(1)));
    double a_post = (1 - x) * a + x * sum(theta.col(1))
      * sum(theta.col(0) / theta.col(1));
    double A_post = x * tau * sum(theta.col(1));
    mu = R::rnorm(a_post, sqrt(A_post));
  }
  if (!fix_t) {
    double w_post = w + theta.n_rows;
    double W_post = W + sum(pow(theta.col(0) - mu, 2.0) / theta.col(1));
    tau = 1.0 / R::rgamma(w_post / 2.0, 2.0 / W_post);
  }
  arma::vec out(2);
  out(0) = mu;
  out(1) = tau;
  return out;
}
