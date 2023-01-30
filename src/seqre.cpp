#include "base.h"
#include "seq_measures.h"
#include "base_measures.h"
#if defined(_OPENMP)
#include <omp.h>
#endif

Rcpp::List marg_cpp(arma::vec y, arma::uword k, arma::uvec z0,
                    seq_measure& seq, base_measure& bas,
                    arma::uword burn, arma::uword thin) {
  arma::uword n = y.n_elem;
  arma::cube recs(k, n, bas.dim + 1);
  arma::mat recs_bas(k, bas.n_pars);
  arma::mat recs_seq(k, seq.n_pars);
  int rec_cnt = 0;
  arma::uvec nums(n, arma::fill::ones);
  nums = cumsum(nums) - 1;
  arma::uvec z = z0;
  arma::uvec z_uq = sort(unique(z));
  arma::uvec cnts = arma::hist(z, nums);

  arma::vec bpars = bas.pars();
  arma::vec bhpars = bas.hpars();
  arma::vec spars = seq.pars();
  arma::vec shpars = seq.hpars();

  arma::mat theta(n, bas.dim);
  for (arma::uword i = 0; i < n; i++) {
    theta.row(i) = bas.rtheta(y(i)).t();
  }
  for (arma::uword iter = 0; iter < burn + thin * k; iter++) {
    arma::uvec z_uq = unique(z);
    arma::uvec cnts = hist(z, nums);
    arma::uword n_uq = z_uq.n_elem;
    arma::mat theta_uq = theta.rows(z_uq);
    bpars = bas.update(theta_uq);

    for (arma::uword i = 0; i < n; i++) {
      cnts(z(i))--;
      if (cnts(z(i)) == 0) {
        n_uq--;
      }
      double q0 = seq.p_new(n - 1, n_uq) * bas.q(y(i));
      arma::vec qj(n, arma::fill::zeros);
      for (arma::uword j = 0; j < n; j++) {
        if (j == i) {
          qj(j) = q0;
        } else {
          qj(j) = seq.p_old(n - 1, n_uq, cnts(z(j))) / cnts(z(j)) *
            bas.q(y(i), theta.row(z(j)).t());
        }
      }
      qj = qj / sum(qj);
      arma::uword z_prop = sample(nums, 1, true, qj)(0);
      if (z_prop == i) {
        theta.row(i) = bas.rtheta(y(i)).t();
        z(i) = z_prop;
        n_uq++;
      } else {
        z(i) = z(z_prop);
      }
      cnts(z(i))++;
    }
    spars = seq.update(n, arma::hist(z));
    if ((iter >= burn) && (iter % thin == 0)) {
      recs.slice(0).row(rec_cnt) = arma::conv_to<arma::vec>::from(z + 1).t();
      recs.slice(1).row(rec_cnt) = theta.col(0).t();
      recs.slice(2).row(rec_cnt) = theta.col(1).t();
      recs_bas.row(rec_cnt) = bas.pars().t();
      recs_seq.row(rec_cnt) = seq.pars().t();
      rec_cnt++;
    }
  }
  return Rcpp::List::create(recs, recs_bas, recs_seq);
}

base_measure* get_bas(arma::uword bas_idx, arma::vec bas_pars,
                      arma::vec bas_hpars) {
  base_measure* out(NULL);
  switch (bas_idx) {
    case 0: {
      g_normls* bas = new g_normls();
      bas->mu = bas_pars(0);
      bas->tau = bas_pars(1);
      bas->s = bas_pars(2);
      bas->S = bas_pars(3);
      bas->a = bas_hpars(0);
      bas->A = bas_hpars(1);
      bas->fix_m = (bool)(bas_hpars(2));
      bas->w = bas_hpars(3);
      bas->W = bas_hpars(4);
      bas->fix_t = (bool)(bas_hpars(5));
      out = bas;
      break;
    }
    default: {
      Rcpp::stop("Unsupported base measure.");
    }
  }
  return out;
}

seq_measure* get_seq(arma::uword seq_idx, arma::vec seq_pars,
                      arma::vec seq_hpars) {
  seq_measure* out(NULL);
  switch (seq_idx) {
    case 0: {
      sq_dirichlet* seq = new sq_dirichlet();
      seq->alpha = seq_pars(0);
      seq->c = seq_hpars(0);
      seq->C = seq_hpars(1);
      seq->fix_a = (bool)(seq_hpars(2));
      out = seq;
      break;
    }
    case 1: {
      sq_pitmanyor* seq = new sq_pitmanyor();
      seq->d = seq_pars(0);
      if (seq_pars(0) < 0) {
        seq->alpha = seq_pars(1) * abs(floor(seq_pars(0)));
      } else {
        seq->alpha = seq_pars(1);
      }
      out = seq;
      break;
    }
    case 2: {
      sq_gnedin0* seq = new sq_gnedin0();
      seq->gamma = seq_pars(0);
      out = seq;
      break;
    }
    default: {
      Rcpp::stop("Unsupported base measure.");
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List marg_cpp(arma::vec y, arma::uword k, arma::uvec z0,
                    arma::uword bas_idx,
                    arma::vec bas_pars, arma::vec bas_hpars,
                    arma::uword seq_idx,
                    arma::vec seq_pars, arma::vec seq_hpars,
                    arma::uword burn, arma::uword thin) {
  base_measure* bas = get_bas(bas_idx, bas_pars, bas_hpars);
  seq_measure* seq = get_seq(seq_idx, seq_pars, seq_hpars);
  Rcpp::List res_marg = marg_cpp(y, k, z0, *seq, *bas, burn, thin);
  return res_marg;
}

// [[Rcpp::export]]
Rcpp::List seqre_cpp(Rcpp::List phi, arma::uword n,
                     arma::uword bas_idx,
                     arma::vec bas_pars, arma::vec bas_hpars,
                     arma::uword seq_idx,
                     arma::vec seq_pars, arma::vec seq_hpars,
                     double eps, arma::uword inc, arma::uword max_it) {
  base_measure* bas = get_bas(bas_idx, bas_pars, bas_hpars);
  seq_measure* seq = get_seq(seq_idx, seq_pars, seq_hpars);
  Rcpp::List res_full(phi.length());
  // TODO: Test parallel with OpenMP
  for (int k = 0; k < phi.length(); k++) {
    arma::mat p = phi[k];
    arma::vec cnts_ = p.col(0) * n;
    arma::uvec cnts = arma::conv_to<arma::uvec>::from(cnts_);
    arma::uvec cnts_new = rnext(*seq, inc, cnts, true);
    arma::uword it = 0;
    arma::uword nn = n;
    double rat0 = (double)(arma::max(cnts)) / nn;
    double rat1 = (double)(arma::max(cnts_new)) / (nn + inc);
    double score = abs(rat1 - rat0);
    nn = nn + inc;
    while((score > eps) & (it < max_it)) {
      cnts = cnts_new;
      cnts_new = rnext(*seq, inc, cnts, true);
      it++;
      rat0 = (double)(arma::max(cnts)) / nn;
      rat1 = (double)(arma::max(cnts_new)) / (nn + inc);
      score = abs(rat1 - rat0);
      nn = nn + inc;
    }
    cnts = cnts_new;
    arma::mat theta = p.cols(arma::span(1, bas->dim));
    arma::uword k_new = cnts.n_elem - p.n_rows;
    if (k_new != 0) {
      arma::mat theta_new(k_new, bas->dim);
      for (arma::uword kk = 0; kk < k_new; kk++) {
        theta_new.row(kk) = bas->rtheta().t();
      }
      theta = join_cols(theta, theta_new);
    }
    cnts_ = arma::conv_to<arma::vec>::from(cnts);
    res_full[k] = join_rows(cnts_ / arma::sum(cnts_), theta);
  }
  return res_full;
}
