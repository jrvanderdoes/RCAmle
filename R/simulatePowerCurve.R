
#' Simulate Power Curve
#'
#' Simulate data, test for change using 3 methods, plot power curve
#'
#' @param Us3 Vector of numerics indicating the first 3 parameters
#' @param U4s Vector of numerics indicating the 4th parameters (i.e. multiple to
#'     examine)
#' @param nSims Integer indicating the number of simulations for each trial
#' @param lowerEst Vector of numerics indicating the lower bound for the 4
#'     parameters. Variance cannot be negative and the third term must also be
#'     above 0.
#' @param upperEst  Vector of numerics indicating the upper bound for the 4
#'     parameters.
#' @param alpha Numeric indicating the significance of interest.
#' @param burnin Integer indicating the number of burnin iteration to use
#' @param k Integer indicating the point at when the process should change from
#'     using beta1 to using beta2
#' @param N Integer indicating how many observations the final data should be.
#' @param errorType (Optional) String indicating the type of error. Options are 'Normal',
#'     'Bernoulli', and 'Exponential'. Default is Normal.
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
#' @param trimAmt (Optional) Numeric indicating amount to trim. Default is 4.
#' @param silent (Optional) Boolean indicating if progress should be output.
#'     Default is FALSE
#'
#' @return List with two terms.
#'     \enumerate{
#'         \item Plot of power curve
#'         \item Data.frame of data used for power curve
#'     }
#' @export
#'
#' @examples
#' pc1 <- simulatePowerCurve(Us3 = c(0, 0.5, 0.5), U4s = c(0,0.5),
#'                   nSims = 5, lowerEst = c(-Inf,0,10^-8,-Inf),
#'                   upperEst = rep(Inf,4), alpha = 0.05,
#'                   burnin = 500, k = 0.5 * 100, N = 100,
#'                   errorType = 'Normal')
simulatePowerCurve <- function(Us3, U4s, nSims, lowerEst, upperEst,
                       alpha, burnin, k, N,
                       errorType = 'Normal',
                       par1HetereoLocation=NA, par1HetereoMult=NA,
                       par2HetereoLocation=NA, par2HetereoMult=NA,
                       trimAmt=4, silent=FALSE){

  nStart <- computeTrim(trimAmt,'Start', N)
  nEnd <- computeTrim(trimAmt,'End',N)

  if(is.na(par1HetereoLocation) && is.na(par1HetereoMult) &&
     is.na(par2HetereoLocation) && is.na(par2HetereoMult)){
    skedast <- 'BothHomoSk'
  } else if(!is.na(par1HetereoLocation) && !is.na(par1HetereoMult)){
    skedast <- 'par1HeteroSk'
  }else if(!is.na(par2HetereoLocation) && !is.na(par2HetereoMult)){
    skedast <- 'par2HeteroSk'
  } else if(!is.na(par1HetereoLocation) && !is.na(par1HetereoMult) &&
            !is.na(par2HetereoLocation) && !is.na(par2HetereoMult)){
    skedast <- 'BothHeteroSk'
  }else{
    stop('Some Error')
  }

  heteroWLS <- ifelse(skedast=='BothHomoSk', FALSE, TRUE)

  powerData <- data.frame('ChangeSize'=NA,
                          'MLE'=NA, 'Vost'=NA, 'WLS'=NA)
  st = Sys.time()
  for(j in 1:length(U4s)){
    U4 <- U4s[j]

    estims <- data.frame("TN"=NA, "WLSVal"=NA,
                         "cutoff_ML"=NA, "cutoff_V"=NA,
                         "TNMax"=NA, "WLSMax"=NA)
    for(i in 1:nSims){

      y <- generateRCA1Data(pars = c(Us3,U4), k = k, burnin = burnin,
                        iterations = N, errorType = errorType,
                        par1HetereoLocation = par1HetereoLocation,
                        par1HetereoMult = par1HetereoMult,
                        par2HetereoLocation = par2HetereoLocation,
                        par2HetereoMult = par2HetereoMult,
                        silent=silent)

      estims[i,] <- compute3Methods(u = c(Us3,U4), y = y,
                                    lower = lowerEst, upper = upperEst,
                                    nStart = nStart, nEnd = nEnd,
                                    returnEstims = FALSE, heteroWLS = heteroWLS)

      oneIterTime <- difftime(Sys.time(), st, units = 'min')/((j-1)*nSims+i)
      numIterRem <- (length(U4s)-j)*nSims + (nSims-i)

      if(!silent) cat(paste0('Completed: ', j, '/',length(U4s), " - ",
                             i,'/',nSims, ' (',round(oneIterTime*numIterRem,3),
                             'mins remaining',')','\n'))
    }

    ## Save what I need
    powerData[j,] <- c(U4,
                       sum(estims$TN>estims$cutoff_ML)/nSims,
                       sum(estims$TN>estims$cutoff_V)/nSims,
                       sum(estims$WLSVal>estims$cutoff_ML)/nSims
    )
  }

  # To remove warnings
  ChangeSize <- MLE <- Vost <- WLS <- NULL

  plot <- ggplot2::ggplot(powerData)+
    ggplot2::geom_point(mapping=ggplot2::aes(x=ChangeSize,y=MLE, col='MLE'),
                        alpha=0.5) +
    ggplot2::geom_line(mapping=ggplot2::aes(x=ChangeSize,y=MLE,
                                            col='MLE', linetype='MLE'),
                       size=1) +
    ggplot2::geom_point(mapping=ggplot2::aes(x=ChangeSize,y=Vost,
                                             col='MLE-Vost'), alpha=0.5) +
    ggplot2::geom_line(mapping=ggplot2::aes(x=ChangeSize,y=Vost,
                                            col='MLE-Vost', linetype='MLE-Vost'),
                       size=1)+
    ggplot2::geom_point(mapping=ggplot2::aes(x=ChangeSize,y=WLS, col='WLS'),
                        alpha=0.5) +
    ggplot2::geom_line(mapping=ggplot2::aes(x=ChangeSize,y=WLS, col='WLS',
                                            linetype='WLS'), size=1)+
    ggplot2::geom_hline(mapping = ggplot2::aes(yintercept=alpha),
                        linetype='dashed', alpha=0.5)+
    ggplot2::geom_vline(mapping = ggplot2::aes(xintercept=Us3[1]),
                        linetype='dotted', alpha=0.5) +
    ggplot2::theme_bw() +
    ggplot2::labs(fill  = "Method", color = "Method") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5),
          axis.title = ggplot2::element_text(size=14),
          axis.title.y = ggplot2::element_blank(),
          axis.text = ggplot2::element_text(size=14),
          legend.text = ggplot2::element_text(size=12),
          legend.title = ggplot2::element_text(size=12),
          legend.position = 'bottom') +
    ggplot2::xlab(bquote(delta)) +
    ggplot2::scale_colour_manual(name = "Method:",
                        labels = c('MLE',"MLE-Vost","WLS"),
                        values = RColorBrewer::brewer.pal(3,'Set1')[c(1,2,3)]) +
    ggplot2::scale_linetype_manual(name = "Method:",
                          labels = c('MLE',"MLE-Vost","WLS"),
                          values = c('dashed',"solid", 'longdash'))


  list(plot, powerData)
}
