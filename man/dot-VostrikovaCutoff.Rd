% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computeCutoff.R
\name{.VostrikovaCutoff}
\alias{.VostrikovaCutoff}
\title{Vostrikova Cutoff}
\usage{
.VostrikovaCutoff(alpha, N, nStart, nEnd)
}
\arguments{
\item{alpha}{Numeric indicating the significance of interest.}

\item{N}{Numeric indicating the length of data}

\item{nStart}{Numeric indicating starting value for trim.}

\item{nEnd}{Numeric indicating ending value for trim.}
}
\value{
Numeric indicating the cutoff value
}
\description{
Cufoff given in Theorem 2.1 with appropriate finite data approx dicussed in
Vostrikova.
}
\examples{
# This is an internal function and will not be viewable. See use in
#     .computeCutoff.
}
