
#' Compute TN
#'
#' This (internal) function computes a form of Tn (form 2.5), technically sqrt(TN).
#'  See for example compute3Methods.
#'
#' @param u Vector of numerics for the 4 parameters of interest, c(b1,s1^2,s2^2,b2).
#'     This is an initial estimate
#' @param y Vector of numerics for the data
#' @param N (Optional) Integer indicating length of y. Added to minimize calls.
#'     Default is length(y).
#' @param lower (Optional) Vector of numerics indicating the lower bound. Default
#'     is c(-Inf,0,10^-8,-Inf). Variance cannot be negative and the third term
#'     must also be above 0.
#' @param upper (Optional) Vector of numerics indicating the upper bound. Default
#'     is the maximal bound, c(Inf,Inf,Inf,Inf).
#' @param nStart (Optional) Integer indicating starting point to check data.
#'     Default results in minimal trimming (4 for params).
#' @param nEnd (Optional) Integer indicating ending point to check data.
#'     Default results in minimal trimming (4 for params).
#' @param returnEstims (Optional) Boolean indicating if the estimated parameters
#'     should be returned. Default is FALSE.
#' @param estimVarLS (Optional) Boolean to indicate if we should use LS to
#'     estimate var2. Default is FALSE.
#' @param var1BaseEstim (Optional) Numeric indicating base var1 estimate. Default
#'     is 0.5.
#' @param var2BaseEstim (Optional) Numeric indicating base var1 estimate. Default
#'     is 0.5. May not be used if estimVarLS.
#' @param maxOptimIters (Optional) Numeric indicating number of times to re-try
#'     optimization. Can sometimes add if get randomly odd result. Default is 2.
#'
#' @return Depends on returnEstims
#'     \itemize{
#'         \item If true: List with element 1 being the WLS estimates, element
#'             2 being the MLE estimates, and element 3 being vector with
#'             sqrt Tn estimate, location of max Tn, and LNN.
#'         \item If false: Vector with sqrt Tn estimate and loc of max Tn.
#'     }
#' @noRd
.computeTN <- function(u, y, N=length(y),
                       lower=c(-Inf,0,10^-8,-Inf), upper=c(Inf,Inf,Inf,Inf),
                       nStart=NA, nEnd=NA, returnEstims = FALSE,
                       estimVarLS=FALSE, var1BaseEstim = 0.5,
                       var2BaseEstim = 0.5, maxOptimIters = 2){
  # Trim Values
  nStart <- ifelse(is.na(nStart),computeTrim(4,'Start'),nStart)
  nEnd <- ifelse(is.na(nEnd),computeTrim(4,'End',N),nEnd)

  ## Setup DFs
  data <- data.frame("k"=nStart:nEnd, "LNk"=NA)
  mleEstim <- wlsEstim <- data.frame("b1"=NA,"v1"=NA,"v2"=NA,"b2"=NA)

  # Get starting estimatates
  tmpBeta <- .WLSestimate(y, N, lower, upper)[1]
  tmpEstims <- c(tmpBeta,var1BaseEstim,
                 ifelse(estimVarLS,
                        .LSVarestimate(list(y),tmpBeta),
                        var2BaseEstim))

  # Compute for no change
  LNN <- .estimateMaxLikelihood(y=y, k=N,
                       estimates=tmpEstims,
                       lower=lower, upper=upper,
                       maxOptimIters=maxOptimIters, returnWithVal=TRUE)[4]

  ## Estimate at all candidate CPs
  wlsEstim <- t(vapply(data$k,.WLSestimate, c(1,2,3,4),
                       y=y, lower=lower, upper=upper))
  tmpEstims <- ifelse(estimVarLS,
                      t(mapply(.LSVarestimate,
                               y=list(y), beta1=wlsEstim[,1],
                               beta2=wlsEstim[,2], k=data$k)),
                      var2BaseEstim)
  mleEstim <- t(mapply(.estimateMaxLikelihood_array, y=list(y), k=data$k,
                       estim1=wlsEstim[,1],
                       estim2=var1BaseEstim,
                       estim3=tmpEstims,
                       estim4=wlsEstim[,4],
                       lower=list(lower), upper=list(upper),
                       maxOptimIters=maxOptimIters))
  data$LNk <- mapply(.L_array, u1 = mleEstim[,1],
                     u2 = mleEstim[,2], u3 = mleEstim[,3],
                     u4 = mleEstim[,4], k = data$k, y = list(y),
                     N=N)

  ## Compute sqrt(TN)
  TN <- sqrt(max(2 * (data$LNk - LNN)))

  # Returns (add nStart to shift results properly)
  if(returnEstims)
    return(
      list(wlsEstim[which.max(data$LNk),],
           mleEstim[which.max(data$LNk),],
           c( TN, which.max(data$LNk)+nStart, LNN)
      )
    )

  return(c(TN, which.max(data$LNk)+nStart))
}
