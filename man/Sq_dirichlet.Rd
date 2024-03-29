% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seq_measure.R
\name{Sq_dirichlet}
\alias{Sq_dirichlet}
\title{Dirichlet Sequence Measure.}
\usage{
Sq_dirichlet(alpha = 1, c = NULL, C = NULL)
}
\arguments{
\item{alpha}{The concentration parameter for the Dirichlet process. Must be
greater than 0.}

\item{c}{The prior primary shape parameter for \code{alpha}.}

\item{C}{The prior secondary shape parameter for \code{alpha}.}
}
\value{
A \code{seq_measure} object for use in the exchangeable sequence
resampling scheme for mixtures.
}
\description{
Dirichlet Sequence Measure.
}
\seealso{
\code{\link[=seq_measure]{seq_measure()}}, \code{\link[=seqre]{seqre()}}
}
