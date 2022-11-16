#ifndef __BASE0__
#define __BASE0__

#include <RcppArmadillo.h>


arma::uvec sample(const arma::uvec& x, const int& size, const bool& replace,
                  const arma::vec& probs);

arma::vec rG(double s, double S, double mu, double tau);

arma::vec rG(double y, double s, double S, double mu, double tau);

double q(double y, double alpha, double s, double S, double mu, double tau);

double q(double y, double d, double a, double s, double S, double mu,
         double tau, unsigned int mk);

double q(double y, double mu, double v);

#endif
