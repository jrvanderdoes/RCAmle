% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computeTN.R
\name{.computeTN}
\alias{.computeTN}
\title{Compute TN}
\usage{
.computeTN(
  u,
  y,
  N = length(y),
  lower = c(-Inf, 0, 10^-8, -Inf),
  upper = c(Inf, Inf, Inf, Inf),
  nStart = NA,
  nEnd = NA,
  returnEstims = FALSE,
  estimVarLS = FALSE,
  var1BaseEstim = 0.5,
  var2BaseEstim = 0.5,
  maxOptimIters = 2
)
}
\arguments{
\item{u}{Vector of numerics for the 4 parameters of interest, c(b1,s1^2,s2^2,b2).
This is an initial estimate}

\item{y}{Vector of numerics for the data}

\item{N}{(Optional) Integer indicating length of y. Added to minimize calls.
Default is length(y).}

\item{lower}{(Optional) Vector of numerics indicating the lower bound. Default
is c(-Inf,0,10^-8,-Inf). Variance cannot be negative and the third term
must also be above 0.}

\item{upper}{(Optional) Vector of numerics indicating the upper bound. Default
is the maximal bound, c(Inf,Inf,Inf,Inf).}

\item{nStart}{(Optional) Integer indicating starting point to check data.
Default results in minimal trimming (4 for params).}

\item{nEnd}{(Optional) Integer indicating ending point to check data.
Default results in minimal trimming (4 for params).}

\item{returnEstims}{(Optional) Boolean indicating if the estimated parameters
should be returned. Default is FALSE.}

\item{estimVarLS}{(Optional) Boolean to indicate if we should use LS to
estimate var2. Default is FALSE.}

\item{var1BaseEstim}{(Optional) Numeric indicating base var1 estimate. Default
is 0.5.}

\item{var2BaseEstim}{(Optional) Numeric indicating base var1 estimate. Default
is 0.5. May not be used if estimVarLS.}

\item{maxOptimIters}{(Optional) Numeric indicating number of times to re-try
optimization. Can sometimes add if get randomly odd result. Default is 2.}
}
\value{
Depends on returnEstims
\itemize{
\item If true: List with element 1 being the WLS estimates, element
2 being the MLE estimates, and element 3 being vector with
sqrt Tn estimate, location of max Tn, and LNN.
\item If false: Vector with sqrt Tn estimate and loc of max Tn.
}
}
\description{
This (internal) function computes a form of Tn (form 2.5), technically sqrt(TN)
}
\examples{
# This is an internal function. See for example compute3Methods.
}
