% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute3Methods.R
\name{compute3Methods}
\alias{compute3Methods}
\title{Compute 3 Methods}
\usage{
compute3Methods(
  u,
  y,
  lower = c(-Inf, 0, 10^-8, -Inf),
  upper = c(Inf, Inf, Inf, Inf),
  alpha = 0.05,
  nStart = NA,
  nEnd = NA,
  returnEstims = FALSE,
  heteroWLS = FALSE
)
}
\arguments{
\item{u}{Vector of numerics for the 4 parameters of interest, c(b1,s1^2,s2^2,b2).
This is an initial estimate}

\item{y}{Vector of numerics for the data}

\item{lower}{(Optional) Vector of numerics indicating the lower bound. Default
is c(-Inf,0,10^-8,-Inf). Variance cannot be negative and the third term
must also be above 0.}

\item{upper}{(Optional) Vector of numerics indicating the upper bound. Default
is the maximal bound, c(Inf,Inf,Inf,Inf).}

\item{alpha}{(Optional) Numeric indicating the significance of interest.
Default is 0.05.}

\item{nStart}{(Optional) Integer indicating starting point to check data.
Default results in minimal trimming (4 for params).}

\item{nEnd}{(Optional) Integer indicating ending point to check data.
Default results in minimal trimming (4 for params).}

\item{returnEstims}{(Optional) Boolean indicating if the estimated parameters
should be returned. Default is FALSE.}

\item{heteroWLS}{(Optional) Boolean indicating if the data is hetereoskedastic.
This changes the WLS estimation. Default is FALSE.}
}
\value{
Depends on returnEstims
\itemize{
\item If true: List with element 1 being the WLS estimates, element
2 being the MLE estimates, and element 3 being the same as if false.
\item If false: Vector with Tn estimate, WLS estimate, cutoff for MLE/WLS,
cutoff for MLE_Vost, loc of max Tn, loc of max WLS.
}
}
\description{
This function computes the three methods (WLS, MLE, MLE_Vost). It has some
optimization to improve speed over calling the three individual methods.
}
\examples{
data <- generateRCA1Data(pars=c(0.5,1,0.25,-0.5), k=250, burnin=1000,
                        iterations=500)
compute3Methods(c(0.4,0.9,0.1,-0.1), data)
}
