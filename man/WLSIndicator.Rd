% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indicators.R
\name{WLSIndicator}
\alias{WLSIndicator}
\title{Weighted Least Squares Estimate Indicator}
\usage{
WLSIndicator(y, N = length(y), alpha = 0.05, heteroWLS = FALSE, trim = log(N))
}
\arguments{
\item{y}{Vector of numerics for the data}

\item{N}{(Optional) Integer indicating length of data. Added for speed
increases.}

\item{alpha}{(Optional) Numeric indicating the significance of interest.
Default is 0.05.}

\item{heteroWLS}{(Optional) Boolean indicating if the data is assumed to be
hetereoskedastic or not.}

\item{trim}{(Optional) Numeric indicating trim amount}
}
\value{
Vector with values estimate, cutoff, and loc of max
}
\description{
This function
}
\examples{
data <- generateRCA1Data(pars=c(0.5,1,0.25,-0.5), k=250, burnin=1000,
                        iterations=500)
WLSIndicator(data)
}
