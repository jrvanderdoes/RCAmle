
#' Estimate Maximum Likelihood
#'
#' An (internal) function that optimizes the functions for parameters based on
#'     the data with the changepoint. See use in computeTN.
#'
#' @param y Vector of numerics for the data
#' @param k Integer indicating the point at when the process should change from
#'     using beta1 to using beta2
#' @param estimates Vector of 4 numerics indicating initial est
#' @param lower Vector of numerics indicating the lower bound. Variance cannot
#'     be negative and the third term must also be above 0.
#' @param upper Vector of numerics indicating the upper bound.
#' @param maxOptimIters (Optional) Numeric indicating number of times to re-try
#'     optimization. Can sometimes add if get randomly odd result. Default is 2.
#' @param returnWithVal (Optional) Boolean indicating if the estimated parameters
#'     should be returned with the estimated values.
#'
#' @return Depends on returnWithVal
#'     \itemize{
#'         \item If TRUE: Vector of the estimated parameters and the estimated
#'             value using the parameters.
#'         \item If FALSE: Vector of the estimated parameters
#'     }
#' @noRd
.estimateMaxLikelihood <- function(y, k, estimates, lower, upper,
                          maxOptimIters=2, returnWithVal=FALSE){

  notDone <- TRUE
  optimIters <- 0

  if(length(estimates)==4){
    ## If it gets an extreme value, try re-running
    while(notDone && optimIters < maxOptimIters){
      ov <- stats::optim(par=estimates, k=k, y=y, N=length(y),
                  fn = .L, method = "L-BFGS-B",#"CG",
                  control = list(fnscale = -1),#'CG'
                  lower = lower,  upper = upper)

      if((ov$par[1]!=lower[1] && ov$par[1]!=upper[1]) &&
         (ov$par[2]!=lower[2] && ov$par[2]!=upper[2]) &&
         (ov$par[3]!=lower[3] && ov$par[3]!=upper[3]) &&
         (ov$par[4]!=lower[4] && ov$par[4]!=upper[4])){
        notDone = FALSE
      }
      optimIters <- optimIters + 1
    }

    # if(optimIters == maxOptimIters && maxOptimIters>1){
    #   warning(paste("Potentially Poor Convergence after",
    #                 maxOptimIters,"iterations"))
    # }
  } else if(length(estimates)==3){

    ## If it gets an extreme value, try re-running
    while(notDone && optimIters < maxOptimIters){
      ov <- stats::optim(par=estimates, y=y,
                  fn = function(u,y){sum(.l(y, c(u[1], u[2], u[3])))},
                  method = "L-BFGS-B",#"CG",
                  control = list(fnscale = -1),
                  lower = lower,  upper = upper)

      if((ov$par[1]!=lower[1] && ov$par[1]!=upper[1]) &&
         (ov$par[2]!=lower[2] && ov$par[2]!=upper[2]) &&
         (ov$par[3]!=lower[3] && ov$par[3]!=upper[3])){
        notDone = FALSE
      }
      optimIters <- optimIters + 1
    }

  }

  if(returnWithVal)
    return(c(ov$par, ov$value))
  return(ov$par)
}


#' Estimate Maximum Likelihood - Array
#'
#' An (internal) function that optimizes the functions for parameters based on
#'     the data with the changepoint. Organizes to use as array for speed. See
#'     use in computeTN.
#'
#' @param y List of vector of numerics for the data
#' @param k Integer indicating the point at when the process should change from
#'     using estim1 to using estim4
#' @param estim1 Vector of numerics for estimates of beta1
#' @param estim2 Vector of numerics for estimates of sigma1^2
#' @param estim3 Vector of numerics for estimates of sigma2^2
#' @param estim4 Vector of numerics for estimates of beta2
#' @param lower List of vector of numerics indicating the lower bound. Variance
#'     cannot be negative and the third term must also be above 0.
#' @param upper List of vector of numerics indicating the upper bound.
#' @param maxOptimIters (Optional) Numeric indicating number of times to re-try
#'     optimization. Can sometimes add if get randomly odd result. Default is 2.
#' @param returnWithVal (Optional) Boolean indicating if the estimated parameters
#'     should be returned with the estimated values.
#'
#' @return Depends on returnWithVal
#'     \itemize{
#'         \item If TRUE: Vector of the estimated parameters and the estimated
#'             value using the parameters.
#'         \item If FALSE: Vector of the estimated parameters
#'     }
#' @noRd
.estimateMaxLikelihood_array <- function(y, k, estim1, estim2, estim3, estim4,
                                         lower, upper, maxOptimIters=2,
                                         returnWithVal=FALSE){
  # Organize inputs
  y <- unlist(y)
  lower <- unlist(lower)
  upper <- unlist(upper)
  estimates <- c(estim1,estim2,estim3,estim4)

  notDone <- TRUE
  optimIters <- 0

  ## If it gets an extreme value, try re-running
  while(notDone && optimIters < maxOptimIters){
    ov <- stats::optim(par=estimates, k=k, y=y, N=length(y),
                fn = .L, method = "L-BFGS-B",#"CG",
                control = list(fnscale = -1),#'CG'
                lower = lower,  upper = upper)

    if((ov$par[1]!=lower[1] && ov$par[1]!=upper[1]) &&
       (ov$par[2]!=lower[2] && ov$par[2]!=upper[2]) &&
       (ov$par[3]!=lower[3] && ov$par[3]!=upper[3]) &&
       (ov$par[4]!=lower[4] && ov$par[4]!=upper[4])){
      notDone = FALSE
    }
    optimIters <- optimIters + 1
  }

  # if(optimIters == maxOptimIters){
  #   warning(paste("Potentially Poor Convergence after",
  #                 maxOptimIters,"iterations"))
  # }

  if(returnWithVal)
    return(c(ov$par, ov$value))
  return(ov$par)
}
