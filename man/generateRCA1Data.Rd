% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateRCA1Data.R
\name{generateRCA1Data}
\alias{generateRCA1Data}
\title{Generate Data}
\usage{
generateRCA1Data(
  pars,
  k,
  burnin,
  iterations,
  errorType = "Normal",
  stationarySims = 10000,
  stationaryCutoff = -0.01,
  par1HetereoLocation = NA,
  par1HetereoMult = NA,
  par2HetereoLocation = NA,
  par2HetereoMult = NA,
  silent = F
)
}
\arguments{
\item{pars}{Vector of parameters for the errors. The format is always similar,
c(beta1, p1, p2, beta2); however, p1 and p2 depend on the errorType.
Errors are always centered around 0
\itemize{
\item Normal: beta1, sigma1^2, sigma2^2, beta2
\item Bernoulli: beta1, prob1, prob2, beta2
\item Exponential: beta1, rate1, rate2, beta2
}}

\item{k}{Integer indicating the point at when the process should change from
using beta1 to using beta2}

\item{burnin}{Integer indicating the number of burn in iterations to use. If the
model is non-stationary, then no burn ins will be used.}

\item{iterations}{Integer indicating how many observations the final data
should be.}

\item{errorType}{String indicating the type of error. Options are 'Normal',
'Bernoulli', and 'Exponential'.}

\item{stationarySims}{Integer indicating the number of simulations to be used
to check for stationarity.}

\item{stationaryCutoff}{Numeric indicating the cutoff for stationarity.
Stationarity is used to determine if burnin should be used. This numeric
can allow for nonstationary regimes, jsut be careful the results don't
blow up.}

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

\item{silent}{(Optional) Boolean indicating if output is requested}
}
\value{
Vector with the sequence of RCA data of length iterations according
to the params
}
\description{
This function generates data according to an RCA(1) process.
}
\examples{
data <- generateRCA1Data(pars=c(0.5,1,0.25,-0.5), k=250, burnin=1000,
                        iterations=500)

data1 <- generateRCA1Data(pars=c(1.1,1,0.25,-0.5), k=250, burnin=1000,
                         iterations=500)

data2 <- generateRCA1Data(pars=c(1.1,1,0.25,-0.5), k=250, burnin=1000,
                         iterations=500,
                         par1HetereoLocation=100,
                         par1HetereoMult=2)
}
