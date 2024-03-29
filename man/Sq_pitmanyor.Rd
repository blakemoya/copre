% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seq_measure.R
\name{Sq_pitmanyor}
\alias{Sq_pitmanyor}
\title{Pitman-Yor Sequence Measure.}
\usage{
Sq_pitmanyor(d, alpha = 1, m = 1L)
}
\arguments{
\item{d}{The discount parameter for the Pitman-Yor process. Must be less than
1.}

\item{alpha}{The concentration parameter for the Pitman-Yor process. Must be
greater than -\code{sigma} if \code{sigma} is in [0, 1), else ignored.}

\item{m}{A positive integer used to set \code{theta = m * abs(sigma)} if \code{sigma}
is negative.}
}
\value{
A \code{seq_measure} object for use in the exchangeable sequence
resampling scheme for mixtures.
}
\description{
Pitman-Yor Sequence Measure.
}
\seealso{
\code{\link[=seq_measure]{seq_measure()}}, \code{\link[=seqre]{seqre()}}
}
