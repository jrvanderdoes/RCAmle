#' Generate Data
#'
#' This function generates data according to an RCA(1) process.
#'
#' @param pars Vector of parameters for the errors. The format is always similar,
#'     c(beta1, p1, p2, beta2); however, p1 and p2 depend on the errorType.
#'     Errors are always centered around 0
#'     \itemize{
#'         \item Normal: beta1, sigma1^2, sigma2^2, beta2
#'         \item Bernoulli: beta1, prob1, prob2, beta2
#'         \item Exponential: beta1, rate1, rate2, beta2
#'     }
#' @param k Integer indicating the point at when the process should change from
#'     using beta1 to using beta2
#' @param burnin Integer indicating the number of burnin iteration to use
#' @param iterations Integer indicating how many observations the final data
#'     should be.
#' @param errorType String indicating the type of error. Options are 'Normal',
#'     'Bernoulli', and 'Exponential'.
#' @param stationarySims Integer indicating the number of simulations to be used
#'     to check for stationarity.
#' @param stationaryCutoff Numeric indicating the cutoff for stationarity.
#'     Stationarity is used to determine if burnin should be used. This numeric
#'     can allow for nonstationary regimes, jsut be careful the results don't
#'     blow up.
#' @param par1HetereoLocation (Optional) This can be used to allow
#'     hetereoskedasticity in the par1 (2nd term in pars) error. This determines
#'     location for it.
#' @param par1HetereoMult (Optional) This can be used to allow
#'     hetereoskedasticity in the par1 (2nd term in pars) error. This determines
#'     the amount to multiply par1 by for it.
#' @param par2HetereoLocation (Optional) This can be used to allow
#'     hetereoskedasticity in the par2 (3rd term in pars) error. This determines
#'     location for it.
#' @param par2HetereoMult (Optional) This can be used to allow
#'     hetereoskedasticity in the par2 (3rd term in pars) error. This determines
#'     the amount to multiply par2 by for it.
#' @param silent (Optional) Boolean indicating if output is requested
#'
#' @return Vector with the sequence of RCA data of length iterations according
#'     to the params
#' @export
#'
#' @examples
#' data <- generateRCA1Data(pars=c(0.5,1,0.25,-0.5), k=500, burnin=1000,
#'                         iterations=1000)
#'
#' data1 <- generateRCA1Data(pars=c(1.1,1,0.25,-0.5), k=500, burnin=1000,
#'                          iterations=1000)
#'
#' data2 <- generateRCA1Data(pars=c(1.1,1,0.25,-0.5), k=500, burnin=1000,
#'                          iterations=1000,
#'                          par1HetereoLocation=0.2,
#'                          par1HetereoMult=2)
generateRCA1Data <- function(pars, k, burnin, iterations,
                         errorType = 'Normal',
                         stationarySims = 10000,
                         stationaryCutoff = -0.01,
                         par1HetereoLocation=NA, par1HetereoMult=NA,
                         par2HetereoLocation=NA, par2HetereoMult=NA,
                         silent=F){
  ## Setup
  par1HetereoLIters <- par1HetereoLocation * iterations
  par2HetereoLIters <- par2HetereoLocation * iterations

  ## Generate Data
  if(errorType == 'Normal'){
    ## Check for stationary
    if(.checkNonStationary(pars, stationarySims, stationaryCutoff)){
      burnin <- 0
      if(!silent) warning("Model nonstationary - No burnin used")
    }

    ## Model Data Generations
    H1 <- .dataGen(pars=pars[c(1,2,3)], burnin=burnin, iterations=k,
                  errorType=errorType,
                  par1HetereoLocation=par1HetereoLIters,
                  par1HetereoMult=par1HetereoMult,
                  par2HetereoLocation=par2HetereoLIters,
                  par2HetereoMult=par2HetereoMult)
    H2 <- .dataGen(pars=pars[c(4,2,3)], burnin=0, iterations=iterations-k,
                  errorType=errorType,
                  par1HetereoLocation=par1HetereoLIters-k,
                  par1HetereoMult=par1HetereoMult,
                  par2HetereoLocation=par2HetereoLIters-k,
                  par2HetereoMult=par2HetereoMult,
                  start = H1[length(H1)])

  }else{
    ## Check for stationary
    if(!silent) warning("Cannot verify stationarity for this error type")

    ## Model Data Generations
    H1 <- .dataGen(pars=pars[c(1,2,3)], burnin=burnin, iterations=k,
                  errorType=errorType,
                  par1HetereoLocation=par1HetereoLIters,
                  par1HetereoMult=par1HetereoMult,
                  par2HetereoLocation=par2HetereoLIters,
                  par2HetereoMult=par2HetereoMult)
    H2 <- .dataGen(pars=pars[c(4,2,3)], burnin=0, iterations=iterations-k,
                  errorType=errorType,
                  par1HetereoLocation=par1HetereoLIters-k,
                  par1HetereoMult=par1HetereoMult,
                  par2HetereoLocation=par2HetereoLIters-k,
                  par2HetereoMult=par2HetereoMult,
                  start = H1[length(H1)])
  }

  # Combine and return two regimes
  c(H1, H2)
}



#' Check Non stationarity
#'
#' This (internal) function check if parameters (for normally distributed errors)
#'     are in the stationary or nonstationary regime.
#'
#' Upcoming: Check when non-normal. Simple mathematically, just need to get
#'     around to it.
#'
#' @param pars Vector of parameters for the normal errors,
#'     c(beta1,sigma1^2,sigma2^2,beta2). Only beta1 and sigma1^2 are used.
#' @param stationarySims Integer for the number of simulations used to check
#'     stationarity.
#' @param stationaryCutoff Numeric indicating the cutoff to determine
#'     stationarity. True cutoff would be 0, but due to finite samples or wanting
#'     to permit some nonstationarity or only want very stationary process, can
#'     make influence user to choose different number.
#'
#' @return Boolean indicating if the process is nonstationary (T) or stationary (F)
#'
#' @examples
#' # This is an internal function and will not be viewable. See use in
#' #     generateData.
.checkNonStationary <- function(pars, stationarySims, stationaryCutoff){
  ## Check for stationary
  #   stationary = E log|beta+e1|<0
  #     E log |pars[1]+N(0,sqrt(pars[2]))|
  stationaryValue <-
    mean(
      log(
        abs(pars[1] + rnorm(stationarySims, mean=0, sd = sqrt(pars[2])))
      )
    )

  stationaryValue > stationaryCutoff
}


#' Data Generation
#'
#' This (internal) function generates data according to an RCA process, using
#'     the given parameters and information.
#'
#' @param pars Vector of parameters for the errors. The format is always similar,
#'     c(beta1, p1, p2); however, p1 and p2 depend on the errorType.
#'     Errors are always centered around 0
#'     \itemize{
#'         \item Normal: beta1, sigma1^2, sigma2^2
#'         \item Bernoulli: beta1, prob1, prob2
#'         \item Exponential: beta1, rate1, rate2
#'     }
#' @param burnin Integer indicating the number of burnin iteration to use
#' @param iterations Integer indicating how many observations the final data
#'     should be.
#' @param errorType String indicating the type of error. Options are 'Normal',
#'     'Bernoulli', and 'Exponential'.
#' @param par1HetereoLocation (Optional) This can be used to allow
#'     hetereoskedasticity in the par1 (2nd term in pars) error. This determines
#'     location for it.
#' @param par1HetereoMult (Optional) This can be used to allow
#'     hetereoskedasticity in the par1 (2nd term in pars) error. This determines
#'     the amount to multiply par1 by for it.
#' @param par2HetereoLocation (Optional) This can be used to allow
#'     hetereoskedasticity in the par2 (3rd term in pars) error. This determines
#'     location for it.
#' @param par2HetereoMult (Optional) This can be used to allow
#'     hetereoskedasticity in the par2 (3rd term in pars) error. This determines
#'     the amount to multiply par2 by for it.
#' @param start (Optional) Numeric, if not NA, will start from this value
#'
#' @return Vector with the sequence of RCA data of length iterations according
#'     to the params
#'
#' @examples
#' # This is an internal function and will not be viewable. See use in
#' #     generateData.
.dataGen  <- function(pars, burnin, iterations, errorType,
                     par1HetereoLocation, par1HetereoMult,
                     par2HetereoLocation, par2HetereoMult,
                     start=NA){

  ## Pars
  b <- pars[1]

  iters <- burnin + iterations
  par1HetereoLoc <- burnin + par1HetereoLocation
  par2HetereoLoc <- burnin + par2HetereoLocation

  ## Build
  if(errorType=='Normal'){
    sigma1 <- sqrt(pars[2])
    sigma2 <- sqrt(pars[3])

    if(!is.na(par1HetereoLoc) &&
       par1HetereoLoc <= iters){
      eps1 <- c(
        rnorm(par1HetereoLoc, mean=0,
              sd=sigma1),
        rnorm(iters-par1HetereoLoc, mean=0,
              sd=sigma1*sqrt(par1HetereoMult))
      )
    } else{

      eps1 <- rnorm(iters, mean=0, sd=sigma1)
    }

    if(!is.na(par2HetereoLoc) &&
       par2HetereoLoc <= iters){
      eps2 <- c(
        rnorm(par2HetereoLoc, mean=0,
              sd=sigma2),
        rnorm(iters-par2HetereoLoc, mean=0,
              sd=sigma2*sqrt(par2HetereoMult))
      )
    }else{
      eps2 <- rnorm(iters, mean=0, sd=sigma2)
    }

  }else if(errorType=='Bernoulli'){

    if(!is.na(par1HetereoLoc) &&
       par1HetereoLoc <= iters){
      eps1 <- c(
        rbinom(n=par1HetereoLoc, size=1,
               prob = pars[2]) - pars[2],
        rbinom(n=iters-par1HetereoLoc, size=1,
               prob = pars[2]*par1HetereoMult) - pars[2]*par1HetereoMult
      )

    } else{
      eps1 <- rbinom(n=iters, size=1, prob = pars[2]) - pars[2]
    }

    if(!is.na(par2HetereoLoc) &&
       par2HetereoLoc <= iters){
      eps2 <- c(
        rbinom(n=par2HetereoLoc, size=1,
               prob = pars[3]) - pars[3],
        rbinom(n=iters-par2HetereoLoc, size=1,
               prob = pars[3]*par2HetereoMult) - pars[3]*par2HetereoMult
      )
    }else{
      eps2 <- rbinom(n=iters, size=1, prob = pars[3]) - pars[3]
    }

  }else if(errorType=='Exponential'){

    if(!is.na(par1HetereoLoc) &&
       par1HetereoLoc <= iters){
      eps1 <- c(
        rexp(n=par1HetereoLoc,
             rate=pars[2]) - 1/pars[2],
        rexp(n=iters-par1HetereoLoc,
             rate=pars[2]*par1HetereoMult) - 1/(pars[2]*par1HetereoMult)
      )
    } else{
      eps1 <- rexp(n=iters,rate=pars[2]) - 1/pars[2]
    }

    if(!is.na(par2HetereoLoc) &&
       par2HetereoLoc <= iters){
      eps2 <- c(
        rexp(n=par2HetereoLoc,
             rate=pars[3]) - 1/pars[3],
        rexp(n=iters-par2HetereoLoc,
             rate=pars[3]*par2HetereoMult) - 1/(pars[3]*par2HetereoMult)
      )
    } else{
      eps2 <- rexp(n=iters,rate=pars[3]) - 1/pars[3]
    }

  }else{
    stop('Error Distribution Incorrectly Specified')
  }

  ## Determine if starting or continuing
  y <- rep(NA,iters)
  if(is.na(start))
    y[1] <- eps2[1]
  else
    y[1] <- (b+eps1[1])*start+eps2[1]

  ## Generate Data
  for(i in 2:iters){
    y[i] <- (b+eps1[i])*y[i-1]+eps2[i]
  }

  ## Drop burnin if needed
  retVal <- NA

  if(burnin>0)
    retVal <- y[-c(1:burnin)]
  else
    retVal <- y

  retVal
}
