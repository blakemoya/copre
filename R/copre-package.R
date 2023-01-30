#' CopRe Tools for Nonparametric Martingale Posterior Sampling
#'
#' Performs Bayesian nonparametric density estimation using Martingale posterior
#' distributions including the Copula Resampling (CopRe) algorithm. Also
#' included are a Gibbs sampler for the marginal Gibbs-type mixture model and an
#' extension to include full uncertainty quantification via a predictive
#' sequence resampling (SeqRe) algorithm. The CopRe and SeqRe samplers generate
#' random nonparametric distributions as output, leading to complete
#' nonparametric inference on posterior summaries. Routines for calculating
#' arbitrary functionals from the sampled distributions are included as well as
#' an important algorithm for finding the number and location of modes, which
#' can then be used to estimate the clusters in the data using, for example,
#' k-means. Implements work developed in Moya B., Walker S. G. (2022).
#'
#' @aliases copre-package
#' @docType package
#' @author Blake Moya <blakemoya@utexas.edu>
#' @references \itemize{
#'  \item Fong, E., Holmes, C., Walker, S. G. (2021). Martingale Posterior
#'  Distributions. arXiv. DOI: \doi{10.48550/arxiv.2103.15671}
#'  \item Moya B., Walker S. G. (2022). Uncertainty Quantification and the
#'  Marginal MDP Model. arXiv. DOI: \doi{10.48550/arxiv.2206.08418}
#'  \item Escobar M. D., West, M. (1995) Bayesian Density Estimation and
#'  Inference Using Mixtures. Journal of the American Statistical Association.
#'  DOI: \doi{10.1080/01621459.1995.10476550}
#' }
#' @import Rcpp pracma abind dirichletprocess
#' @importFrom Rcpp sourceCpp
#' @importFrom pracma gradient
#' @importFrom abind abind
#' @importFrom utils packageVersion
#' @useDynLib copre
#' @keywords internals
"_PACKAGE"
