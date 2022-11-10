// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// chk_omp
bool chk_omp();
RcppExport SEXP _copre_chk_omp() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(chk_omp());
    return rcpp_result_gen;
END_RCPP
}
// copre_cpp
arma::mat copre_cpp(arma::mat y_perm, arma::vec alpha, double rho, arma::uword N, arma::uword k, arma::mat& P, arma::vec grd, arma::uword nthreads);
RcppExport SEXP _copre_copre_cpp(SEXP y_permSEXP, SEXP alphaSEXP, SEXP rhoSEXP, SEXP NSEXP, SEXP kSEXP, SEXP PSEXP, SEXP grdSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y_perm(y_permSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type grd(grdSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(copre_cpp(y_perm, alpha, rho, N, k, P, grd, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// mdp_cpp
Rcpp::List mdp_cpp(arma::vec y, arma::uword k, double alpha, double mu, double tau, double s, double S, double c, double C, double a, double A, double w, double W, bool fix_a, bool fix_m, bool fix_t, arma::uword burn, arma::uword thin);
RcppExport SEXP _copre_mdp_cpp(SEXP ySEXP, SEXP kSEXP, SEXP alphaSEXP, SEXP muSEXP, SEXP tauSEXP, SEXP sSEXP, SEXP SSEXP, SEXP cSEXP, SEXP CSEXP, SEXP aSEXP, SEXP ASEXP, SEXP wSEXP, SEXP WSEXP, SEXP fix_aSEXP, SEXP fix_mSEXP, SEXP fix_tSEXP, SEXP burnSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::uword >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type C(CSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type W(WSEXP);
    Rcpp::traits::input_parameter< bool >::type fix_a(fix_aSEXP);
    Rcpp::traits::input_parameter< bool >::type fix_m(fix_mSEXP);
    Rcpp::traits::input_parameter< bool >::type fix_t(fix_tSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(mdp_cpp(y, k, alpha, mu, tau, s, S, c, C, a, A, w, W, fix_a, fix_m, fix_t, burn, thin));
    return rcpp_result_gen;
END_RCPP
}
// polyaurn_cpp
Rcpp::List polyaurn_cpp(arma::cube theta, arma::mat hyper, double eps, double ups, arma::uword nthreads);
RcppExport SEXP _copre_polyaurn_cpp(SEXP thetaSEXP, SEXP hyperSEXP, SEXP epsSEXP, SEXP upsSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type hyper(hyperSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double >::type ups(upsSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(polyaurn_cpp(theta, hyper, eps, ups, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// mpy_cpp
Rcpp::List mpy_cpp(arma::vec y, arma::uword k, double d, double alpha, double mu, double tau, double s, double S, double c, double C, double a, double A, double w, double W, bool fix_a, bool fix_m, bool fix_t, arma::uword burn, arma::uword thin);
RcppExport SEXP _copre_mpy_cpp(SEXP ySEXP, SEXP kSEXP, SEXP dSEXP, SEXP alphaSEXP, SEXP muSEXP, SEXP tauSEXP, SEXP sSEXP, SEXP SSEXP, SEXP cSEXP, SEXP CSEXP, SEXP aSEXP, SEXP ASEXP, SEXP wSEXP, SEXP WSEXP, SEXP fix_aSEXP, SEXP fix_mSEXP, SEXP fix_tSEXP, SEXP burnSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::uword >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type C(CSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type W(WSEXP);
    Rcpp::traits::input_parameter< bool >::type fix_a(fix_aSEXP);
    Rcpp::traits::input_parameter< bool >::type fix_m(fix_mSEXP);
    Rcpp::traits::input_parameter< bool >::type fix_t(fix_tSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(mpy_cpp(y, k, d, alpha, mu, tau, s, S, c, C, a, A, w, W, fix_a, fix_m, fix_t, burn, thin));
    return rcpp_result_gen;
END_RCPP
}
// polyaurngen_cpp
Rcpp::List polyaurngen_cpp(arma::umat theta, arma::vec alpha, double eps, double ups, arma::uword nthreads);
RcppExport SEXP _copre_polyaurngen_cpp(SEXP thetaSEXP, SEXP alphaSEXP, SEXP epsSEXP, SEXP upsSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::umat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double >::type ups(upsSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(polyaurngen_cpp(theta, alpha, eps, ups, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// evalmdpolya_cpp
arma::mat evalmdpolya_cpp(Rcpp::List& mdpolya, arma::vec& x, arma::uword f, arma::uword nthreads);
RcppExport SEXP _copre_evalmdpolya_cpp(SEXP mdpolyaSEXP, SEXP xSEXP, SEXP fSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type mdpolya(mdpolyaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type f(fSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(evalmdpolya_cpp(mdpolya, x, f, nthreads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_copre_chk_omp", (DL_FUNC) &_copre_chk_omp, 0},
    {"_copre_copre_cpp", (DL_FUNC) &_copre_copre_cpp, 8},
    {"_copre_mdp_cpp", (DL_FUNC) &_copre_mdp_cpp, 18},
    {"_copre_polyaurn_cpp", (DL_FUNC) &_copre_polyaurn_cpp, 5},
    {"_copre_mpy_cpp", (DL_FUNC) &_copre_mpy_cpp, 19},
    {"_copre_polyaurngen_cpp", (DL_FUNC) &_copre_polyaurngen_cpp, 5},
    {"_copre_evalmdpolya_cpp", (DL_FUNC) &_copre_evalmdpolya_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_copre(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
