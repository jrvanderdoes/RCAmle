#' Compute 3 Methods
#'
#' This function computes the three methods (WLS, MLE, MLE_Vost). It has some
#'     optimization to improve speed over calling the three individual methods.
#'
#' @param u Vector of numerics for the 4 parameters of interest, c(b1,s1^2,s2^2,b2).
#'     This is an initial estimate
#' @param y Vector of numerics for the data
#' @param lower (Optional) Vector of numerics indicating the lower bound. Default
#'     is c(-Inf,0,10^-8,-Inf). Variance cannot be negative and the third term
#'     must also be above 0.
#' @param upper (Optional) Vector of numerics indicating the upper bound. Default
#'     is the maximal bound, c(Inf,Inf,Inf,Inf).
#' @param alpha (Optional) Numeric indicating the significance of interest.
#'     Default is 0.05.
#' @param nStart (Optional) Integer indicating starting point to check data.
#'     Default results in minimal trimming (4 for params).
#' @param nEnd (Optional) Integer indicating ending point to check data.
#'     Default results in minimal trimming (4 for params).
#' @param returnEstims (Optional) Boolean indicating if the estimated parameters
#'     should be returned. Default is FALSE.
#' @param heteroWLS (Optional) Boolean indicating if the data is hetereoskedastic.
#'     This changes the WLS estimation. Default is FALSE.
#'
#' @return Depends on returnEstims
#'     \itemize{
#'         \item If true: List with element 1 being the WLS estimates, element
#'             2 being the MLE estimates, and element 3 being the same as if false.
#'         \item If false: Vector with Tn estimate, WLS estimate, cutoff for MLE/WLS,
#'             cutoff for MLE_Vost, loc of max Tn, loc of max WLS.
#'     }
#' @export
#'
#' @examples
#' data <- generateRCA1Data(pars=c(0.5,1,0.25,-0.5), k=500, burnin=1000,
#'                         iterations=1000)
#' compute3Methods(c(0.4,0.9,0.1,-0.1), data)
compute3Methods <- function(u, y,
                            lower=c(-Inf,0,10^-8,-Inf), upper=c(Inf,Inf,Inf,Inf),
                            alpha = 0.05, nStart = NA, nEnd=NA,
                            returnEstims = FALSE, heteroWLS = FALSE){
  # Params
  N <- length(y)

  # Compute Tn Value (For MLE orig and vost)
  results <- .computeTN(u=u, y=y, N=N,
                       lower=lower, upper=upper,
                       nStart=nStart, nEnd=nEnd,
                       returnEstims = returnEstims)

  # Compute WLS
  if(is.na(nStart)){
    trim <- trimSelect(4,'Start')-1
  }else{
    trim <- nStart-1
  }

  if(heteroWLS)
    PreVal <- computeWLSVals_Hetero(y=y, N=N, trim=trim)
  else
    PreVal <- computeWLSVals(y=y, N=N, trim=trim)

  # Get Cutoffs
  cutoff_ML <- .computeCutoff(alpha, N, "General") # MLE/WLS Method
  cutoff_V <- .computeCutoff(alpha=alpha, N=N, "Vostrikova") # Vostrikova Method

  # Return Values
  if(returnEstims)
    return(
      list("WLS"=results[[1]],
           "MLE"=results[[2]],
           c(results[[3]][1],PreVal[1],cutoff_ML,
             cutoff_V,results[[3]][2],PreVal[2])
      )
    )

  return(c(results[1],PreVal[1],cutoff_ML,
           cutoff_V,results[2],PreVal[2]))
}
