.onAttach <- function(libname, pkgname) {
  packageStartupMessage('CopRe v', utils::packageVersion(pkgname),
                        ': OpenMP ', ifelse(chk_omp(), 'enabled', 'disabled'))
  register_autoplot_s3_methods()
}
