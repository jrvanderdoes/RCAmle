
#' Compute Weighted Least Square Values
#'
#' This uses weighted least squares to estimate paramters for RCA(1) with changes.
#'     This is for homoskedastic data. This work is replication of that done
#'     in Horvath, Trapani 2022.
#'
#' @param y Vector of numerics for the data
#' @param N (Optional) Numeric indicating length of the data. Included to
#'     increase speed. Default is length(y).
#' @param trim (Optional) Numeric indicating trim amount. Default is log(N).
#'
#' @return Vector of two values: value and location of max
#' @export
#'
#' @examples
#' data <- generateRCA1Data(pars=c(0.5,1,0.25,-0.5), k=500, burnin=1000,
#'                         iterations=1000)
#' computeWLSVals(data)
computeWLSVals <- function(y, N=length(y), trim=log(N)){

  # Prep and Drop y's (to use lags)
  ylag    <- c(NA,y[-N])
  ylag2   <- ylag^2
  ylagden <- 1+ylag2

  # Weight Function
  w05 <- (((1:N)/N)*(1-(1:N)/N))^0.5

  #  Create num and den for WLS estimation
  num     <- (y*ylag)/ylagden
  den     <- ylag2/ylagden
  deneta  <- ylag2/ylagden^2

  # Estimate eta
  betahat   <- sum(num[-1])/sum(den[-1])
  residuals <- y - betahat*ylag
  res2      <- residuals^2
  a1        <- sum(res2[-1]*deneta[-1])/(N-1)
  a2        <- sum(den[-1])/(N-1)
  a2        <- a2^2
  ## Page 10 - 3.6
  etasq     <- a1/a2
  eta       <- sqrt(etasq)

  beta1 <- beta2 <- qt <- rep(0,N)

  for(k0 in 2:N){
    num1      <- sum(num[2:k0])
    den1      <- sum(den[2:k0])
    beta1[k0] <- num1/den1
  }

  for(k0 in 1:(N-1)){
    num2      <- sum(num[(k0+1):N])
    den2      <- sum(den[(k0+1):N])
    beta2[k0] <- num2/den2
  }

  for(k0 in 1:N){
    qt[k0] <- sqrt(N)*(k0/N)*(1-k0/N)*(beta1[k0]-beta2[k0])

    if (k0 <= ceiling(trim+1) ||
        k0 > N-ceiling(trim+1)
    ) # +1 for array indexing
      qt[k0] <- 0

  }

  #  Weighted Functionals
  q05 <- abs(qt)/w05
  #  sup norms
  mq05 <- max(q05[-N])/eta

  return(c(mq05, which.max(q05[-N])))
}


#' Compute Weighted Least Square Values - Hetereoskedasticity
#'
#' This uses weighted least squares to estimate paramters for RCA(1) with changes.
#'     This is for hetereoskedastic data. This work is replication of that done
#'     in Horvath, Trapani 2022.
#'
#' @param y Vector of numerics for the data
#' @param N (Optional) Numeric indicating length of the data. Included to
#'     increase speed. Default is length(y).
#' @param trim (Optional) Numeric indicating trim amount. Default is log(N).
#'
#' @return Vector of two values: value and location of max
#' @export
#'
#' @examples
#' data <- generateRCA1Data(pars=c(1.1,1,0.25,-0.5), k=500, burnin=1000,
#'                          iterations=1000,
#'                          par1HetereoLocation=100,
#'                          par1HetereoMult=2)
#' computeWLSVals(data)
computeWLSVals_Hetero <- function(y, N=length(y), trim=log(N)){

  # Prep and Drop y's (to use lags)
  ylag    = c(NA,y[-N])
  ylag2   = ylag^2
  ylagden = 1+ylag2

  # Weight Function
  w05 <- (((1:N)/N)*(1-(1:N)/N))^0.5

  #  Create num and den for WLS estimation
  num     = (y*ylag)/ylagden;
  den     = ylag2/ylagden;
  deneta  = ylag2/ylagden^2;

  # Estimate eta
  betahat   <- sum(num[-1])/sum(den[-1])
  residuals <- y - betahat*ylag
  denb <- (residuals*ylag/ylagden)^2
  bt <- rep(0,N)
  c1t <- rep(0,N)

  for(tk in 2:N){
    bt[tk] <- sum(denb[2:tk])/N
    c1t[tk]<- sum(den[2:tk])/N
  }

  c2t <- c1t[N]-c1t
  c1 <- c1t[N]
  b1 <- bt[N]
  b11<- bt[N-1]
  g0t <- (c1^2)*bt-2*c1*(c1t*bt) + b1*(c1t*c1t)
  g0t <- ifelse(abs(g0t)<10^-16,0,g0t)
  #erdos1 <- c1t*c2t
  #erdos2 <- g0t*erdos1
  erdos <- sqrt(g0t)

  beta1 <- rep(0,N)
  beta2 <- rep(0,N)
  qt <- rep(0,N)

  for(k0 in 2:N){
    num1      = sum(num[2:k0])
    den1      = sum(den[2:k0])
    beta1[k0] = num1/den1
  }

  for(k0 in 1:(N-1)){
    num2      = sum(num[(k0+1):N])
    den2      = sum(den[(k0+1):N])
    beta2[k0] = num2/den2
  }

  for(k0 in 2:N){
    qt[k0] = sqrt(N)*c1t[k0]*c2t[k0]*(beta1[k0]-beta2[k0])

    if (k0 <= ceiling(trim+1) ||
        k0 > N-ceiling(trim+1)
    ) # +1 for array indexing
      qt[k0] = 0
  }

  ## Page 15
  #  Weighted Functionals
  q05 = abs(qt)/erdos
  #  sup norms
  useQ05 <- q05[-c(1,N)]
  dropIdx <- 1
  if(length(which(is.na(useQ05)))){
    dropIdx <- which(is.na(useQ05))
    useQ05 <- useQ05[-dropIdx]
    dropIdx <- dropIdx + 1
  }
  mq05 = max(useQ05)

  return(c(mq05, which.max(useQ05)+dropIdx))
}


#' Weighted Least Squares Estimate
#'
#' @param y Vector of numerics for the data
#' @param k Integer(s) location of change point. Can do this as a vector
#' @param lower Vector of numerics indicating the lower bound. Variance cannot
#'     be negative and the third term must also be above 0.
#' @param upper Vector of numerics indicating the upper bound.
#'
#' @return WLS estimates for beta1, var1, var2, beta2 at potentially different k
#'
#' @examples
#' # Internal function. See for use example binarySegmentationCPDetection.
.WLSestimate <- function(y,k, lower, upper){

  means <- .WLSMeansEstimate(y,k)
  vars <- .WLSVarsEstim(y, k, means)

  result <- c(means[1],vars[1], vars[2], means[2])

  return(
    mapply(function(x,upper, lower){
      if(is.na(x)) return(x)
      if(x>upper) x <- upper
      if(x<lower) x <- lower
      x
    }, x=result, upper=upper, lower=lower)
  )

}


#' Weighted Least Squares Mean Estimates
#'
#' Compute weighted least squares mean estimates.
#'
#' @param y Vector of numerics for the data
#' @param k Integer location of change point
#'
#' @return Vector of two numerics indicating beta1 and beta2 estimates
#'
#' @examples
#' # Internal function, see use in .WLSestimate.
.WLSMeansEstimate <- function(y,k){

  N <- length(y)

  ## Solve For Beta1 and Beta2
  # Regular Case
  if(k!=N){

    denom1 <- 0
    num1 <- 0
    for(i in 2:k){
      denom1 <- denom1 + y[i-1]^2/(1+y[i-1]^2)
      num1 <- num1 + y[i]*y[i-1]/(1+y[i-1]^2)
    }

    denom2 <- 0
    num2 <- 0
    for(i in (k+1):N){
      denom2 <- denom2 + y[i-1]^2/(1+y[i-1]^2)
      num2 <- num2 + y[i]*y[i-1]/(1+y[i-1]^2)
    }


    b1Hat <- num1/denom1
    b2Hat <- num2/denom2

  }else{ # For LNN
    denom1 <- 0
    num1 <- 0
    for(i in 2:k){
      denom1 <- denom1 + y[i-1]^2/(1+y[i-1]^2)
      num1 <- num1 + y[i]*y[i-1]/(1+y[i-1]^2)
    }

    b1Hat <- num1/denom1
    b2Hat <- NA

  }

  return(c(b1Hat,b2Hat))
}


#' Weighted Least Squares Variance Estimates
#'
#' Compute weighted least squares variance estimates. See for example books like
#'     Random Coefficient Autoregressive Models: An Introduction (Des F. Nicholls,
#'     Barry G. Quinn) or A New Iterative Procedure for Estimation of RCA
#'     Parameters Based on Estimating Functions (Norlli Anida Abdullah et al.)
#'
#' @param y List of vector of numerics of data
#' @param changepoint Integer location of change point
#' @param means Vector of two numerics indicating means. Likely from
#'     .WLSMeansEstimate.
#'
#' @return Vector of two numerics indicating var1 and var2 estimates
#'
#' @examples
#' # Internal function, see use in .WLSestimate.
.WLSVarsEstim <- function(y, changepoint, means){

  # Helpful Defns
  n <- length(y)

  if(changepoint < n){
    y1 <- y[1:changepoint]
    y2 <- y[(changepoint+1):n]

    n1 <- length(y1)
    n2 <- length(y2)
  }


  # Resids
  if(changepoint < n){
    u1 <- y1[-1] - means[1] * y1[-n1]
    u2 <- y2[-1] - means[2] * y2[-n2]

    u <- c(u1,u2)

    zBar1 <- sum(y1[-n1]^2)/(n1-1)
    zBar2 <- sum(y2[-n2]^2)/(n2-1)

    zBar <- c(rep(zBar1,length(u1)),rep(zBar2,length(u2)))
  }else{
    u <- y[-1] - means[1] * y[-n]

    zBar <- sum(y[-n]^2)/(n-1)
  }

  #zBar1.1 <- sum(y1[-n1]^2)/(n1)
  #zBar2.1 <- sum(y2[-n2]^2)/(n2)
  #zBar.1 <- c(rep(zBar1.1,length(u1)),rep(zBar2.1,length(u2)))

  ## I tested all sorts of shifting to match the paper
  #       11 for var1 and nearly 100 for var2
  #       Most differences tiny. Not worth effort
  #       Using basic as most simple to explain. Best commented out.
  # Var1
  if(changepoint < n){
    var1 <- sum(u^2*(y[-c(1,n)]^2-zBar^2)) / sum(y[-c(1,n)]^4-zBar^2)
    #var1 <- sum(c(NA,u^2)*(y[-c(n)]^2-zBar.1^2), na.rm = TRUE) / sum(y[-c(n)]^4-c(NA,zBar.1^2), na.rm = TRUE)
  }else{
    var1 <- sum(u^2*(y[-n]^2-zBar^2)) / sum(y[-n]^4-zBar^2)
  }

  # Var2
  var2 <- sum(u^2 - var1 * zBar) / (n-1)
  #var2 <- sum(u^2 - var1.2 * zBar) / (n-1)

  # Require Nonnegativity
  #   Found it best to require at this point and not immediately following calculation
  if(is.na(var1) || var1<0)
    var1 <- 0
  if(is.na(var2) || var2<0)
    var2 <- 0

  c(var1, var2)
}


#' Least Squares Variance1 Estimate
#'
#' This (internal) function computes the least squares estimate of variance1
#'
#' @param y List of vector of data
#' @param beta1 Numeric for beta1 estimate
#' @param beta2 (Optional) Numeric for beta2 estimate
#' @param k (Optional) Integer location of change point when beta2 is present
#'
#' @return Numeric indicating the var1 estimate
#'
#' @examples
#' # This is an internal function. See for example binarySegmentationCPDetection.
.LSVarestimate <- function(y, beta1, beta2=NA, k=NA){

  y <- unlist(y)
  useY <- y[-1]

  if(is.na(beta2))
    resultVal <- sum((useY - beta1 * y[-length(y)])^2)
  else
    resultVal <- sum((useY[1:(k-2)]  - beta1 * y[1:(k-2)])^2) +
    sum((useY[(k-1):length(useY)] - beta2 * y[(k-1):length(useY)])^2)

  resultVal
}
