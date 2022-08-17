#if defined(_OPENMP)
  #include <omp.h>
#endif

//[[Rcpp::export]]
bool chk_omp() {
  #if defined(_OPENMP)
    return true;
  #else
    return false;
  #endif
}
