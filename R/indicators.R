
#' Maximum Liklihood Estimate Indicator
#'
#' Function to  get MLE estimates for data with the related cutoff.
#'
#' @param u Vector of numerics for the 4 parameters of interest, c(b1,s1^2,s2^2,b2).
#'     This is an initial estimate
#' @param y Vector of numerics for the data
#' @param lower Vector of numerics indicating the lower bound.Variance cannot be
#'     negative and the third term must also be above 0.
#' @param upper Vector of numerics indicating the upper bound.
#' @param alpha (Optional) Numeric indicating the significance of interest.
#'     Default is 0.05.
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
#'             sqrt Tn estimate, cutoff, location of max Tn, and LNN.
#'         \item If false: Vector with sqrt Tn estimate, cutoff, and loc of max Tn.
#'     }
#' @export
#'
#' @examples
#' data <- generateRCA1Data(pars=c(0.5,1,0.25,-0.5), k=500, burnin=1000,
#'                         iterations=1000)
#' MLEIndicator(c(0.4,0.9,0.1,-0.1), data, c(-Inf,0,10^-8,-Inf),
#'              c(Inf,Inf,Inf,Inf))
MLEIndicator <- function(u, y,
                         lower, upper,
                         alpha = 0.05,
                         nStart = NA, nEnd=NA,
                         returnEstims = FALSE,
                         estimVarLS=FALSE, var1BaseEstim = 0.5,
                         var2BaseEstim = 0.5,
                         maxOptimIters = 2){

  N <- length(y)

  results <- .computeTN(u=u, y=y, N=N,
                       lower=lower, upper=upper,
                       nStart=nStart, nEnd=nEnd,
                       returnEstims = returnEstims,
                       estimVarLS=estimVarLS, var1BaseEstim = var1BaseEstim,
                       var2BaseEstim = var2BaseEstim,
                       maxOptimIters = maxOptimIters)

  cutoff <- .computeCutoff(alpha, N, "General")

  if(returnEstims)
    return(
      list(results[[1]],
           results[[2]],
           c( results[[3]][1], cutoff, results[[3]][2],results[[3]][3] )
      )
    )

  return(c( results[1], cutoff, results[2] ))
}


#' Maximum Liklihood Estimate with Finite Sample Approximation Indicator
#'
#' Function to  get MLE estimates for data with the related cutoff based on
#'     finite sample approximation.
#'
#' @param u Vector of numerics for the 4 parameters of interest, c(b1,s1^2,s2^2,b2).
#'     This is an initial estimate
#' @param y Vector of numerics for the data
#' @param lower Vector of numerics indicating the lower bound.Variance cannot be
#'     negative and the third term must also be above 0.
#' @param upper Vector of numerics indicating the upper bound.
#' @param alpha (Optional) Numeric indicating the significance of interest.
#'     Default is 0.05.
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
#'             sqrt Tn estimate, cutoff, location of max Tn, and LNN.
#'         \item If false: Vector with sqrt Tn estimate, cutoff, and loc of max Tn.
#'     }
#' @export
#'
#' @examples
#' data <- generateRCA1Data(pars=c(0.5,1,0.25,-0.5), k=500, burnin=1000,
#'                         iterations=1000)
#' MLEVostIndicator(c(0.4,0.9,0.1,-0.1), data, c(-Inf,0,10^-8,-Inf),
#'              c(Inf,Inf,Inf,Inf))
MLEVostIndicator <- function(u, y,
                                lower, upper,
                                alpha = 0.05,
                                nStart = NA, nEnd=NA,
                                returnEstims=FALSE,
                                estimVarLS=FALSE, var1BaseEstim = 0.5,
                                var2BaseEstim = 0.5,
                                maxOptimIters = 2){

  N <- length(y)

  results <- .computeTN(u=u, y=y, N=N,
                       lower=lower, upper=upper,
                       nStart=nStart, nEnd=nEnd,
                       returnEstims = returnEstims,
                       estimVarLS=estimVarLS, var1BaseEstim = var1BaseEstim,
                       var2BaseEstim = var2BaseEstim,
                       maxOptimIters = maxOptimIters)

  U_Star <- ifelse(returnEstims, results[[3]][1], results[1])

  cutoff <- .computeCutoff(alpha, N,"Vostrikova")

  if(returnEstims){
    return(list(results[[1]],
                results[[2]],
                c( U_Star, cutoff, results[[3]][2])))
  }
  return(c( U_Star, cutoff, results[2]))
}


#' Weighted Least Squares Estimate Indicator
#'
#' This function
#'
#' @param y Vector of numerics for the data
#' @param N (Optional) Integer indicating length of data. Added for speed
#'     increases.
#' @param alpha (Optional) Numeric indicating the significance of interest.
#'     Default is 0.05.
#' @param heteroWLS (Optional) Boolean indicating if the data is assumed to be
#'     hetereoskedastic or not.
#' @param trim (Optional) Numeric indicating trim amount
#'
#' @return Vector with values estimate, cutoff, and loc of max
#' @export
#'
#' @examples
#' data <- generateRCA1Data(pars=c(0.5,1,0.25,-0.5), k=500, burnin=1000,
#'                         iterations=1000)
#' WLSIndicator(data)
WLSIndicator <- function(y, N=length(y), alpha=0.05,
                         heteroWLS = FALSE, trim = log(N)){

  #  Darling Erdos Critical Values
  crvalue <- getCutoff(alpha=alpha, N=N, type="General")

  if(heteroWLS)
    retVals <- computeWLSVals_Hetero(y=y, N=N, trim = trim)
  else
    retVals <- computeWLSVals(y=y, N=N, trim = trim)

  c(retVals[1], crvalue, retVals[2])
}
