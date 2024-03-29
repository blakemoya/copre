% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modes.R
\name{antimodes}
\alias{antimodes}
\title{Antiode Extractor}
\usage{
antimodes(obj, mean = FALSE, grd = NULL, idx = FALSE)
}
\arguments{
\item{obj}{A \code{copre_result} or \code{mdp_result} object.}

\item{mean}{A logical value indicating whether to extract the modes of the
mean density of each of the individual sampled density.}

\item{grd}{For \code{mdpolya_result}, a grid on which to evaluate the object.}

\item{idx}{A logical value indicating whether to also return the index within
\code{grd} of the discovered modes.}
}
\value{
A matrix of antimodes values in the support of the \code{copre_result}
density
}
\description{
Extracts the antimodes from a \code{copre_result} or \code{seqre_result}
object.
}
