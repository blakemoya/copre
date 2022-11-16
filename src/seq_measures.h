#ifndef __SEQ__
#define __SEQ__

#include "base.h"

arma::uvec setdiff(arma::uvec x, const arma::uvec& y);

class seq_measure {
public:
  arma::uword n_pars;
  arma::uword n_hpars;
  seq_measure(arma::uword n_pars, arma::uword n_hpars) :
    n_pars(n_pars), n_hpars(n_hpars) { }

  virtual bool validate() = 0;
  virtual double p_new(arma::uword n, arma::uword k) = 0;
  virtual double p_old(arma::uword n, arma::uword k, arma::uword kj) = 0;
  virtual arma::vec pars() = 0;
  virtual arma::vec hpars() = 0;
  virtual arma::vec update(arma::uword n, arma::uvec cnts) = 0;
};

arma::uvec rnext(seq_measure& s, arma::uword n, arma::uvec z_prev, bool tab);

arma::uvec rnext(seq_measure& s, arma::uword n, bool tab);

class sq_gnedin0 : public seq_measure {
public:
  sq_gnedin0() : seq_measure(1, 0) { }

  double gamma = 0;

  bool validate();
  double p_new(arma::uword n, arma::uword k);
  double p_old(arma::uword n, arma::uword k, arma::uword kj);
  arma::vec pars();
  arma::vec hpars();
  arma::vec update(arma::uword n, arma::uvec cnts);
};

class sq_pitmanyor : public seq_measure {
public:
  sq_pitmanyor() : seq_measure(2, 0) { }

  double d = 0;
  arma::uword m = 1;
  double alpha = 1;

  bool validate();
  double p_new(arma::uword n, arma::uword k);
  double p_old(arma::uword n, arma::uword k, arma::uword kj);
  arma::vec pars();
  arma::vec hpars();
  arma::vec update(arma::uword n, arma::uvec cnts);
};

class sq_dirichlet : public seq_measure {
public:
  sq_dirichlet() : seq_measure(1, 2) { }

  double alpha = 1;

  double c = 1;
  double C = 1;
  bool fix_a = false;

  bool validate();
  double p_new(arma::uword n, arma::uword k);
  double p_old(arma::uword n, arma::uword k, arma::uword kj);
  arma::vec pars();
  arma::vec hpars();
  arma::vec update(arma::uword n, arma::uvec cnts);
};

#endif
