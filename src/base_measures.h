#ifndef __BASE__
#define __BASE__

#include <RcppArmadillo.h>

class base_measure {
public:
  arma::uword dim;
  arma::uword n_pars;
  arma::uword n_hpars;
  base_measure(arma::uword dim, arma::uword n_pars, arma::uword n_hpars) :
    dim(dim), n_pars(n_pars), n_hpars(n_hpars) { }

  virtual bool validate() = 0;
  virtual arma::vec rtheta() = 0;
  virtual arma::vec rtheta(double y) = 0;
  virtual double q(double y) = 0;
  virtual double q(double y, arma::vec theta) = 0;
  virtual arma::vec pars() = 0;
  virtual arma::vec hpars() = 0;
  virtual arma::vec update(arma::mat theta) = 0;
};

class g_normls : public base_measure {
public:
  g_normls() : base_measure(2, 4, 4) { }

  double mu = 0;
  double tau = 1;
  double s = 1;
  double S = 1;

  double a = 0;
  double A = 1;
  bool fix_m = false;
  double w = 1;
  double W = 1;
  bool fix_t = false;

  bool validate();
  arma::vec rtheta();
  arma::vec rtheta(double y);
  double q(double y);
  double q(double y, arma::vec theta);
  arma::vec pars();
  arma::vec hpars();
  arma::vec update(arma::mat theta);
};

#endif
