#' Copula Resampling
#'
#' @description A function that samples predictive distributions for univariate
#'   continuous data using the bivariate Gaussian copula.
#'
#' @param data The data from which to sample predictive distributions.
#' @param N The number of unobserved data points to resample for each chain.
#' @param k The number of predictive distributions to sample.
#' @param rho A scalar concentration parameter.
#' @param grd_res The number of points on which to evaluate the predictive
#'   distribution.
#' @param nthreads The number of threads to call for parallel execution.
#' @param gpu A logical value indicating whether or not to use the CUDA
#'   implementation of the algorithm.
#' @param gpu_path The path to the CUDA implementation source code.
#' @param gpu_odir A directory to output the compiled CUDA code.
#' @param gpu_seed A seed for the CUDA random variates.
#'
#' @return A `copre_result` object, whose underlying structure is a list which
#'   contains the following components:
#' @references Fong, E., Holmes, C., Walker, S. G. (2021). Martingale Posterior
#'   Distributions. arXiv. DOI: \doi{10.48550/arxiv.2103.15671}
#' @examples
#' res_cop <- copre(rnorm(50), 10, 10, nthreads = 1)
#' @export
copre <- function(data, N, k, rho = 0.91, grd_res = 1000,
                  nthreads = parallel::detectCores(), gpu = FALSE,
                  gpu_path = NULL, gpu_odir = NULL,
                  gpu_seed = 1234) {
  ls <- list(rho, N)
  lens <- sapply(ls, length)
  argmat <- matrix(nrow = max(lens), ncol = 2)
  argmat[, 1] <- rho
  argmat[, 2] <- N
  core <- function(argvec) {
    rho <- argvec[1]
    N <- argvec[2]
    args <- list(data = data, rho = rho, N = N)
    data <- scale(data)
    n <- length(data)
    perms <- matrix(data, nrow = n, ncol = k)
    perms <- apply(perms, 2, sample)
    alpha <- (2 - 1 / 1:(n + N)) * 1 / (1:(n + N) + 1)
    r_x <- range(data)
    rr_x <- diff(r_x) / 10
    grd <- seq(r_x[1] - rr_x, r_x[2] + rr_x, length = grd_res)
    P <- seq(0.1, 0.9, length = grd_res)
    out <- matrix(P, ncol = k, nrow = grd_res)
    if (gpu) {
      rst <- NULL
      if (!is.loaded("copre_gpu")) {
        if (!file.exists(paste(gpu_odir, 'copre_gpu.o', sep = '/'))) {
          if (!dir.exists(gpu_odir)) {
            dir.create(gpu_odir)
          }
          message('Calling `nvcc`...\n')
          nvcc_out <- system2('nvcc',
                              args = paste(gpu_path,
                                           paste('--shared',
                                                 '-odir', gpu_odir,
                                                 '-o copre_gpu.o')),
                              stdout = TRUE)
          file.rename(paste0('copre_gpu', c('.o', '.exp', '.lib')),
                      paste0(gpu_odir, '/copre_gpu', c('.o', '.exp', '.lib')))
          message(paste0('nvcc> ', nvcc_out, '\n'), sep = '')
        }
        message(paste('Loading `copre_gpu.o` from', gpu_odir, '...\n'))
        dyn.load(paste(gpu_odir, 'copre_gpu.o', sep = '/'))
        message('`copre_gpu.o` loaded.\n')
      }
      env <- environment()
      eval(parse(text = '
        t <- system.time(rst <- .C("copre_gpu",
                                    as.numeric(perms),
                                    as.numeric(alpha),
                                    P = as.numeric(out),
                                    as.numeric(grd),
                                    as.numeric(rho),
                                    as.integer(n),
                                    as.integer(N),
                                    as.integer(k),
                                    as.integer(grd_res),
                                    as.integer(gpu_seed)))
        '), envir = env)
      out <- matrix(rst$P, nrow = k, byrow = TRUE)
    } else {
      t <- system.time(out <- copre_cpp(perms, alpha, rho, N, k, out, grd,
                                        nthreads))
    }
    grd <- grd * attr(data, 'scaled:scale') + attr(data, 'scaled:center')
    attr(out, 'func') <- 'distribution'
    attr(out, 'grid') <- grd
    attr(out, 'args') <- args
    attr(out, 'time') <- t
    class(out) <- c('copre_result', 'grideval_result')
    return(out)
  }
  if (nrow(argmat) == 1) {
    return(core(argmat[1, ]))
  } else {
    return(apply(argmat, 1, core, simplify = FALSE))
  }
}
