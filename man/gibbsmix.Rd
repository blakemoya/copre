% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gibbsmix.R
\name{gibbsmix}
\alias{gibbsmix}
\title{Marginal Gibbs-type Mixture Model Sampler}
\usage{
gibbsmix(data, k, b_msr, s_msr, burn = 1000, thin = 150)
}
\arguments{
\item{data}{The data from which to sample predictive distributions.}

\item{k}{The number of predictive samples to draw.}

\item{b_msr}{A \code{base_measure} object.}

\item{s_msr}{A \code{seq_measure} object.}

\item{burn}{The number of initial sampling iterations to discard, will be
truncated if a non-integer.}

\item{thin}{The number of sampling iterations to discard between records,
will be truncated if a non-integer.}
}
\value{
A \code{seqre_result} object.
}
\description{
A function that samples marginal mixture densities via a
marginal Gibbs sampler.
}
\seealso{
\code{\link[=seqre]{seqre()}}, \code{\link[=seq_measure]{seq_measure()}}, \code{\link[=base_measure]{base_measure()}}
}
