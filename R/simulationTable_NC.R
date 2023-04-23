
#' Simulation Table with Null Case
#'
#' This function is used to generate a table of null rejection frequencies
#'
#' @param betas Vector of numerics indicating betas to examine (remember this is
#'     the null case so only one beta used each iteration)
#' @param varProbRates Vector of 2 values indicating the variance or probabilitise
#'     used in parameters, c(u2, u3).
#' @param nSims Integer indicating the number of simulations for each trial
#' @param iterations Integer indicating how many observations the final data
#'     should be.
#' @param burnin Integer indicating the number of burnin iteration to use
#' @param lowerEstim Vector of numerics indicating the lower bound for the 4
#'     parameters. Variance cannot be negative and the third term must also be
#'     above 0.
#' @param upperEstim Vector of numerics indicating the upper bound for the 4
#'     parameters.
#' @param alpha Numeric indicating the significance of interest.
#' @param errorTypes Vector indicating the errorTypes to investigate. Options
#'     are 'Normal', 'Bernoulli', and 'Exponential'.
#' @param CPLoc Numeric indicating change point location, between 0 and 1.
#' @param seed (Optional) Integer to indicate seed to set at each interation.
#'     Default of NA sets no seed.
#' @param silent (Optional) Boolean indicating if progress should be output.
#'     Default is FALSE
#'
#' @return Data.frame of null rejection frequencies
#' @export
#'
#' @examples
#' data_er <- simulationTable_NC(betas = c(0.5),
#'                               varProbRates = c(0.5, 0.5),
#'                               nSims=5, iterations = c(50),
#'                               burnin = 1000, lowerEst=c(-Inf,0,10^-8,-Inf),
#'                               upperEst=c(Inf,Inf,Inf,Inf), alpha = 0.05,
#'                               errorTypes = c('Normal'),
#'                               CPLoc=0.5)
#'
#' \dontrun{
#' data_er <- simulationTable_NC(betas = c(0.5,0.75,1,1.05),
#'                               varProbRates = c(0.5, 0.5),
#'                               nSims=500, iterations = c(100, 200, 400, 800),
#'                               burnin = 1000, lowerEst=c(-Inf,0,10^-8,-Inf),
#'                               upperEst=c(Inf,Inf,Inf,Inf), alpha = 0.05,
#'                               errorTypes = c('Normal','Bernoulli','Exponential'),
#'                               CPLoc=0.5, seed=1234)
#' }
simulationTable_NC <- function(betas, varProbRates, nSims, iterations, burnin,
                               lowerEstim, upperEstim, alpha, errorTypes,
                               CPLoc=0.5, seed=NA, silent=F,
                               trimAmt = 4){

  sol <- expand.grid(
    'type'=c('MLE','Vost','WLS'),
    'iters'=iterations,
    'beta'=betas,
    'errType'=errorTypes
  )
  sol$Estim <- NA
  sol$Lower <- NA
  sol$Upper <- NA

  ## No Change Point
  for(et in 1:length(errorTypes)){
    errorType <- errorTypes[et]
    if(!silent) cat(paste0('ErrorType: ',errorType, ' (',et,'/',length(errorTypes),')\n'))

    for(b in 1:length(betas)){
      beta <- betas[b]
      if(!silent) cat(paste0('- Starting Beta: ',beta,' (',b,'/',length(betas),')\n'))

      for(i in 1:length(iterations)){
        iters <- iterations[i]
        if(!silent) cat(paste0('-- Running Iters: ',iters,' (',i,'/',length(iterations),')\n'))

        ## Setup Dataframes
        wlsEstim <- mleEstim <- data.frame("b1"=NA,"v1"=NA,"v2"=NA,"b2"=NA)
        valEstim <- data.frame("TN"=NA, "PreVal"=NA,
                               "cutoff_ML"=NA, "cutoff_V"=NA,
                               "TNMax"=NA, "PreValMax"=NA)

        ## Set seed here so possible to continue if trouble arises
        if(!is.na(seed)) set.seed(seed)

        nStart <- computeTrim(trimAmt,'Start')
        nEnd <- computeTrim(trimAmt,'End',iters)

        for(j in 1:nSims){

          y <- generateRCA1Data(pars=c(beta, varProbRates, beta), k=iters*CPLoc,
                            burnin = burnin, iterations = iters,
                            errorType = errorType)

          returnVal <- compute3Methods(u = rep(NA, 4),
                                           y = y,
                                           lower = lowerEstim,
                                           upper = upperEstim,
                                           nStart = nStart,
                                           nEnd = nEnd,
                                           alpha = alpha,
                                           returnEstims = TRUE,
                                           heteroWLS = FALSE)
          wlsEstim[j,] <- returnVal[[1]]
          mleEstim[j,] <- returnVal[[2]]
          valEstim[j,] <- returnVal[[3]]
        }

        ## Record, Save Data
        data <- list(wlsEstim, mleEstim, valEstim)

        MLEEst <- sum(data[[3]]$TN>data[[3]]$cutoff_ML)/nSims
        VostEst<- sum(data[[3]]$TN>data[[3]]$cutoff_V)/nSims
        WLSEst<- sum(data[[3]]$PreVal>data[[3]]$cutoff_ML)/nSims

        MLESigmaHat <- sqrt(stats::var(data[[3]]$TN>data[[3]]$cutoff_ML))
        VostSigmaHat<- sqrt(stats::var(data[[3]]$TN>data[[3]]$cutoff_V))
        WLSSigmaHat<- sqrt(stats::var(data[[3]]$PreVal>data[[3]]$cutoff_ML))

        sol[sol$type=='MLE' &
              sol$iters == iters &
              sol$beta == beta &
              sol$errType == errorType,c('Estim','Lower','Upper')] <-
          c(round(MLEEst,4),
            max(0,round(MLEEst+stats::qnorm(alpha/2)*MLESigmaHat/sqrt(nSims),4)),
            min(1,round(MLEEst-stats::qnorm(alpha/2)*MLESigmaHat/sqrt(nSims),4))
          )

        sol[sol$type=='Vost' &
              sol$iters == iters &
              sol$beta == beta &
              sol$errType == errorType,c('Estim','Lower','Upper')] <-
          c(round(VostEst,4),
            max(0,round(VostEst+stats::qnorm(alpha/2)*VostSigmaHat/sqrt(nSims),4)),
            min(1,round(VostEst-stats::qnorm(alpha/2)*VostSigmaHat/sqrt(nSims),4))
          )

        sol[sol$type=='WLS' &
              sol$iters == iters &
              sol$beta == beta &
              sol$errType == errorType,c('Estim','Lower','Upper')] <-
          c(round(WLSEst,4),
            max(0,round(WLSEst+stats::qnorm(alpha/2)*WLSSigmaHat/sqrt(nSims),4)),
            min(1,round(WLSEst-stats::qnorm(alpha/2)*WLSSigmaHat/sqrt(nSims),4))
          )

        #saveRDS(data,paste0(savePath,'/Sk-',skedast,'_beta-',beta,
        #                    '_iters-',iters,'_',
        #                    substr(errorType,1,1),".RDS"))
      }
    }
  }

  sol
}
