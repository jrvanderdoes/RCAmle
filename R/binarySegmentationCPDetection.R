
#' Binary Segmentation for Change Point Detection
#'
#' This function uses binary segmentation to detect change points.
#'
#' @param fullData Data.frame with two columns. The first has date and the
#'     second has values
#' @param method String indicating method to use for change point detection.
#'     Options are 'MLE', 'Vostrikova' , and 'WLS'.
#' @param lower Vector of numerics indicating the lower bound. Variance cannot
#'     be negative and the third term must also be above 0.
#' @param upper Vector of numerics indicating the upper bound.
#' @param alpha Numeric indicating the significance of interest.
#' @param estimVarLS (Optional) Boolean to indicate if we should use LS to
#'     estimate var2. Default is FALSE.
#' @param var1BaseEstim (Optional) Numeric indicating base var1 estimate. Default
#'     is 0.5.
#' @param var2BaseEstim (Optional) Numeric indicating base var1 estimate. Default
#'     is NA. May not be used if estimVarLS.
#' @param maxOptimIters (Optional) Numeric indicating number of times to re-try
#'     optimization. Can sometimes add if get randomly odd result. Default is 1.
#' @param trimAmt (Optional) Number indicating the amount to trim. Default is 4.
#' @param verifyCPs (Optional) Boolean indicating if change points should be
#'     rechecked at the end. Default is TRUE.
#' @param silent (Optional) Boolean indicating if output should be supressed.
#'     Default is FALSE
#'
#' @return Data.frame with columns 'stDate' (start date of region), 'enDate'
#'     (end date of region), 'WLS_b1' (beta1 estimate using WLS), 'WLS_v1' (var1
#'     estimate using WLS), 'WLS_v2' (var2 estimate using WLS), 'WLS_b2' (beta2
#'     estimate using WLS), 'MLE_b1' (beta1 estimate using MLE), 'MLE_v1' (var1
#'     estimate using MLE), 'MLE_v2' (var2 estimate using MLE), and 'MLE_b2'
#'     (beta2 estimate using MLE).
#' @export
#'
#' @examples
#' data <- generateRCA1Data(pars=c(0.5,1,0.25,-0.5), k=100, burnin=1000,
#'                         iterations=200)
#' data <- data.frame(1:length(data),data)
#' binarySegmentationCPDetection(data,'Vostrikova',
#'                               lower=c(-Inf,0,10^-8,-Inf),
#'                               upper=c(Inf,Inf,Inf,Inf), alpha=0.05)
binarySegmentationCPDetection <- function(fullData, method,
                                  lower, upper,
                                  alpha,
                                  estimVarLS=FALSE,
                                  var1BaseEstim = 0.5,
                                  var2BaseEstim = NA,
                                  maxOptimIters = 1,
                                  trimAmt = 4,
                                  verifyCPs = TRUE,
                                  silent = FALSE){

  # Prepare
  fullData <- fullData[order(fullData[,1]),]
  onlyData <- fullData[,2]

  # Get change points
  if(!silent) cat('\n--Find CPs --\n')
  CPsVals <- .detectChangePoints(data=onlyData, method=method,
                                lower=lower, upper=upper,
                                alpha=alpha,
                                estimVarLS=estimVarLS,
                                var1BaseEstim = var1BaseEstim,
                                var2BaseEstim = var2BaseEstim,
                                maxOptimIters = maxOptimIters,
                                trimAmt=trimAmt,
                                silent=silent)

  CPsVals <- c(0, CPsVals, length(onlyData))

  # Verify as desired
  if(verifyCPs){
    if(!silent) cat('\n-- Verify Step --\n')
    CPsVals <- .verifyChangePoints(CPsVals=CPsVals, data=onlyData,
                      method=method,
                      lower=lower, upper=upper,
                      alpha=alpha,
                      estimVarLS=estimVarLS,
                      var1BaseEstim = var1BaseEstim,
                      var2BaseEstim = var2BaseEstim,
                      maxOptimIters = maxOptimIters,
                      trimAmt=trimAmt,
                      silent=silent)
  }

  CPs <- data.frame('st'=c(CPsVals[-length(CPsVals)]+1),
                    'en'=c(CPsVals[-1]))


  # Get Parameter Estimates
  wls <- mle <- data.frame("b1"=NA, "v1"=NA, "v2"=NA, "b2"=NA)

  for(i in 1:length(CPs$st)){
    # Add bounds to WLS
    wls[i,] <- .WLSestimate(y=onlyData[CPs$st[i]:CPs$en[i]],
                           k=length(onlyData[CPs$st[i]:CPs$en[i]]),
                           lower=lower, upper=upper)

    if(is.na(wls[i,4])){
      tmpEstim <- ifelse(estimVarLS,
                         .LSVarestimate(
                           list(onlyData[CPs$st[i]:CPs$en[i]]),
                           wls[i,1]),
                         ifelse(is.na(var2BaseEstim),
                                stats::median(
                                  abs(onlyData[CPs$st[i]:CPs$en[i]] -
                                        stats::median(onlyData[CPs$st[i]:CPs$en[i]])))/0.6745,
                                var2BaseEstim)
      )
      mle[i,] <- c(.estimateMaxLikelihood(y=onlyData[CPs$st[i]:CPs$en[i]],
                                 k=length(onlyData[CPs$st[i]:CPs$en[i]]),
                                 estimates=c(wls[i,1],
                                             var1BaseEstim,
                                             tmpEstim),
                                 lower=lower, upper=upper),NA)

    }else{
      stop('This should not happen!!')
    }
  }

  returnDF <- data.frame('stDate'=fullData[CPs$st, 1],
                         'enDate'=fullData[CPs$en, 1],
                         'WLS_b1'= wls$b1, 'WLS_v1'= wls$v1,
                         'WLS_v2'= wls$v2, 'WLS_b2'= wls$b2,
                         'MLE_b1'= mle$b1, 'MLE_v1'= mle$v1,
                         'MLE_v2'= mle$v2, 'MLE_b2'= mle$b2 )

  return(returnDF)
}



#' Detect Change Points
#'
#' See use in binarySegmentationCPDetection.
#'
#' @param data Vector of numerics with values, pre-ordered
#' @param method String indicating method to use for change point detection.
#'     Options are 'MLE', 'Vostrikova' , and 'WLS'.
#' @param lower Vector of numerics indicating the lower bound. Variance cannot
#'     be negative and the third term must also be above 0.
#' @param upper Vector of numerics indicating the upper bound.
#' @param alpha Numeric indicating the significance of interest.
#' @param addAmt (Optional) Numeric to add to CPs so can use this function
#'     recursively.
#' @param estimVarLS (Optional) Boolean to indicate if we should use LS to
#'     estimate var2. Default is FALSE.
#' @param var1BaseEstim (Optional) Numeric indicating base var1 estimate. Default
#'     is 0.5.
#' @param var2BaseEstim (Optional) Numeric indicating base var1 estimate. Default
#'     is 0.5. If NA will compute based on robust statistics. May not be used if
#'     estimVarLS.
#' @param maxOptimIters (Optional) Numeric indicating number of times to re-try
#'     optimization. Can sometimes add if get randomly odd result. Default is 1.
#' @param silent (Optional) Boolean indicating if output should be supressed.
#'     Default is FALSE
#'
#' @return Vector of integers indicating the change points
#' @noRd
.detectChangePoints <-  function(data, method,
                                lower, upper,
                                alpha,
                                addAmt = 0,
                                estimVarLS=FALSE, var1BaseEstim = 0.5,
                                var2BaseEstim = 0.5,
                                maxOptimIters = 1,
                                trimAmt = 4,
                                silent = F){
  # Prepare
  sendVar2 <- var2BaseEstim
  if(is.na(var2BaseEstim)) var2BaseEstim <-
      stats::median( abs(data - stats::median(data)))/0.6745

  nStart <- computeTrim(trimAmt,'Start',length(data))
  nEnd <- computeTrim(trimAmt,'End',length(data))

  if(nStart> nEnd) return()

  # Determine Method
  if(method == 'MLE'){
    result <-
      MLEIndicator(u= rep(NA,4), data,
                   lower=lower, upper=upper,
                   alpha = alpha,
                   nStart = nStart, nEnd=nEnd,
                   returnEstims = FALSE,
                   estimVarLS=estimVarLS,
                   var1BaseEstim = var1BaseEstim,
                   var2BaseEstim = var2BaseEstim,
                   maxOptimIters = maxOptimIters)

  } else if(method == 'Vostrikova'){
    result <-
      MLEVostIndicator(u= rep(NA,4), data,
                          lower=lower, upper=upper,
                          alpha = alpha,
                          nStart = nStart, nEnd=nEnd,
                          returnEstims = FALSE,
                          estimVarLS=estimVarLS,
                          var1BaseEstim = var1BaseEstim,
                          var2BaseEstim = var2BaseEstim,
                          maxOptimIters = maxOptimIters)
  } else if(method == 'WLS'){
    result <-
      WLSIndicator(y=data, N=length(data),
                   alpha = alpha)
  } else {
    stop('Incorrect Method specified')
  }

  # No Change Point Detected
  if(is.na(result[1]) ||
     is.na(result[2]) ||
     result[1] <= result[2]){
    return()
  }

  if(!silent) cat(paste0('ChangePoint Detected (',1+addAmt,'-' ,addAmt+length(data),' at ',
                         addAmt+result[3],'): Segment Data and Re-Search\n'))

  return(c(
    .detectChangePoints(data=data[1:result[3]],
                       method=method,
                       lower=lower, upper=upper,
                       alpha=alpha,
                       addAmt = addAmt,
                       estimVarLS=estimVarLS,
                       var1BaseEstim = var1BaseEstim,
                       var2BaseEstim = sendVar2,
                       maxOptimIters = maxOptimIters,
                       trimAmt=trimAmt,
                       silent=silent),
    result[3],
    .detectChangePoints(data=data[(result[3]+1):length(data)],
                       method=method,
                       lower=lower, upper=upper,
                       alpha=alpha,
                       addAmt = addAmt+result[3],
                       estimVarLS=estimVarLS,
                       var1BaseEstim = var1BaseEstim,
                       var2BaseEstim = sendVar2,
                       maxOptimIters = maxOptimIters,
                       trimAmt=trimAmt,
                       silent=silent)+result[3]))

}

#' Verify Change Points
#'
#' See use in binarySegmentationCPDetection.
#'
#' @param CPsVals Vector of integers indicating the detected change points
#' @param data Vector of numerics indicating data values
#' @param method String indicating method to use for change point detection.
#'     Options are 'MLE', 'Vostrikova' , and 'WLS'.
#' @param lower Vector of numerics indicating the lower bound. Variance cannot
#'     be negative and the third term must also be above 0.
#' @param upper Vector of numerics indicating the upper bound.
#' @param alpha Numeric indicating the significance of interest.
#' @param estimVarLS Boolean to indicate if we should use LS to estimate var2.
#' @param var1BaseEstim Numeric indicating base var1 estimate.
#' @param var2BaseEstim Numeric indicating base var1 estimate.
#' @param maxOptimIters Numeric indicating number of times to re-try
#'     optimization. Can sometimes add if get randomly odd result.
#' @param silent (Optional) Boolean indicating if output should be supressed.
#'     Default is FALSE
#'
#' @return Vector of integers indicating the change points
#' @noRd
.verifyChangePoints <- function(CPsVals, data,
                   method,
                   lower, upper,
                   alpha,
                   estimVarLS,
                   var1BaseEstim,
                   var2BaseEstim,
                   maxOptimIters,
                   trimAmt,
                   silent=F){

  newCPVals <- c()
  var2BaseEstim1 <- var2BaseEstim
  # For no CP, i.e. c(0,length(data))
  if(length(CPsVals)==2) CPsVals <- c(CPsVals,length(data))

  for(i in 1:(length(CPsVals)-2)){
    nStart <- computeTrim(trimAmt,'Start',length(data[(CPsVals[i]+1):CPsVals[i+2]]))
    nEnd <- computeTrim(trimAmt,'End',length(data[(CPsVals[i]+1):CPsVals[i+2]]))

    if(nStart > nEnd) next

    if(is.na(var2BaseEstim))
      var2BaseEstim1 <- stats::median( abs(data[(CPsVals[i]+1):CPsVals[i+2]] -
                                             stats::median(data[(CPsVals[i]+1):CPsVals[i+2]])))/0.6745

    # Determine Method
    if(method == 'MLE'){
      result <-
        MLEIndicator(u= rep(NA,4), y = data[(CPsVals[i]+1):CPsVals[i+2]],
                     lower=lower, upper=upper,
                     alpha = alpha,
                     nStart = nStart, nEnd=nEnd,
                     returnEstims = FALSE,
                     estimVarLS=estimVarLS,
                     var1BaseEstim = var1BaseEstim,
                     var2BaseEstim = var2BaseEstim1,
                     maxOptimIters = maxOptimIters)

    } else if(method == 'Vostrikova'){
      result <-
        MLEVostIndicator(u= rep(NA,4), y = data[(CPsVals[i]+1):CPsVals[i+2]],
                         lower=lower, upper=upper,
                         alpha = alpha,
                         nStart = nStart, nEnd=nEnd,
                         returnEstims = FALSE,
                         estimVarLS=estimVarLS,
                         var1BaseEstim = var1BaseEstim,
                         var2BaseEstim = var2BaseEstim1,
                         maxOptimIters = maxOptimIters)
    } else if(method == 'WLS'){
      result <-
        WLSIndicator(y=data[(CPsVals[i]+1):CPsVals[i+2]],
                     N=length(data[(CPsVals[i]+1):CPsVals[i+2]]),
                     alpha = alpha)
    } else {
      stop('Incorrect Method specified')
    }

    if(!(is.na(result[1]) || is.na(result[2]) || result[1] <= result[2])){
      if(!silent) {
        cat(paste0('ChangePoint Detected (',1+CPsVals[i],'-' ,
                   CPsVals[i]+length(data[(CPsVals[i]+1):CPsVals[i+2]]),
                   ' at ',CPsVals[i]+result[3], '): Segment Data and Re-Search\n'))
      }

      newCPVals <- c(newCPVals, result[3]+CPsVals[i])
    }
  }

  c(0,newCPVals,length(data))
}
