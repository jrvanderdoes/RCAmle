% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulatePowerCurve.R
\name{simulatePowerCurve}
\alias{simulatePowerCurve}
\title{Simulate Power Curve}
\usage{
simulatePowerCurve(
  Us3,
  U4s,
  nSims,
  lowerEst,
  upperEst,
  alpha,
  burnin,
  k,
  N,
  errorType = "Normal",
  par1HetereoLocation = NA,
  par1HetereoMult = NA,
  par2HetereoLocation = NA,
  par2HetereoMult = NA,
  trimAmt = 4,
  silent = FALSE
)
}
\arguments{
\item{Us3}{Vector of numerics indicating the first 3 parameters}

\item{U4s}{Vector of numerics indicating the 4th parameters (i.e. multiple to
examine)}

\item{nSims}{Integer indicating the number of simulations for each trial}

\item{lowerEst}{Vector of numerics indicating the lower bound for the 4
parameters. Variance cannot be negative and the third term must also be
above 0.}

\item{upperEst}{Vector of numerics indicating the upper bound for the 4
parameters.}

\item{alpha}{Numeric indicating the significance of interest.}

\item{burnin}{Integer indicating the number of burnin iteration to use}

\item{k}{Integer indicating the point at when the process should change from
using beta1 to using beta2}

\item{N}{Integer indicating how many observations the final data should be.}

\item{errorType}{(Optional) String indicating the type of error. Options are 'Normal',
'Bernoulli', and 'Exponential'. Default is Normal.}

\item{par1HetereoLocation}{(Optional) This can be used to allow
hetereoskedasticity in the par1 (2nd term in pars) error. This determines
location for it.}

\item{par1HetereoMult}{(Optional) This can be used to allow
hetereoskedasticity in the par1 (2nd term in pars) error. This determines
the amount to multiply par1 by for it.}

\item{par2HetereoLocation}{(Optional) This can be used to allow
hetereoskedasticity in the par2 (3rd term in pars) error. This determines
location for it.}

\item{par2HetereoMult}{(Optional) This can be used to allow
hetereoskedasticity in the par2 (3rd term in pars) error. This determines
the amount to multiply par2 by for it.}

\item{trimAmt}{(Optional) Numeric indicating amount to trim. Default is 4.}

\item{silent}{(Optional) Boolean indicating if progress should be output.
Default is FALSE}
}
\value{
List with two terms.
\enumerate{
\item Plot of power curve
\item Data.frame of data used for power curve
}
}
\description{
Simulate data, test for change using 3 methods, plot power curve
}
\examples{
pc1 <- simulatePowerCurve(Us3 = c(0, 0.5, 0.5), U4s = c(0,0.5),
                  nSims = 5, lowerEst = c(-Inf,0,10^-8,-Inf),
                  upperEst = rep(Inf,4), alpha = 0.05,
                  burnin = 500, k = 0.5 * 100, N = 100,
                  errorType = 'Normal')
}
