#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>


arma::uvec sample(const arma::uvec& x, const int& size, const bool& replace,
                  const arma::vec& probs){
  return Rcpp::RcppArmadillo::sample(x, size, replace, probs);
}

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

double q(double y, double d, double a, double s, double S, double mu,
         double tau, unsigned int mk) {
  double q0 = (a + d * mk) * tgamma((1.0 + s) / 2.0) / (tgamma(s / 2.0) * sqrt(s));
  q0 *= pow(1.0 + pow(y - mu, 2.0) / ((1 + tau) * S), -(1 + s) / 2.0);
  q0 /= sqrt((1 + tau) * S / s);
  return q0;
}

double q(double y, double mu, double v) {
  return exp(-pow(y - mu, 2.0) / (2 * v)) / sqrt(2 * v);
}

