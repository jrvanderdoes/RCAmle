
#' Compute Cutoff
#'
#' An (internal) function to compute the cutoff value for the statistic of
#'  interest. See use in compute3Methods.
#'
#' @param alpha Numeric indicating the significance of interest.
#' @param N Numeric indicating the length of data
#' @param type String indicating type of cutoff. Options are 'General' or
#'     'Vostrikova'.
#' @param nStart (Optional) Numeric indicating starting value for trim. Default
#'     is NA.
#' @param nEnd (Optional) Numeric indicating ending value for trim. Default
#'     is NA.
#'
#' @return Numeric indicating the cutoff value
#' @noRd
.computeCutoff <- function(alpha, N, type, nStart=NA, nEnd=NA){

  nStart <- ifelse(is.na(nStart),computeTrim(4,'Start'),nStart)
  nEnd <- ifelse(is.na(nEnd),computeTrim(4,'End',N),nEnd)

  if(type=="General"){
    cutoff <- .GeneralCutoff(alpha, N)
  }else if(type == "Vostrikova"){
    cutoff <- .VostrikovaCutoff(alpha, N, nStart, nEnd)
  } else{
    stop('Error: Review value of type')
  }

  cutoff
}


#' General Cutoff
#'
#' Cutoff given in Theorem 2.1. See use in .computeCutoff.
#'
#' @param alpha Numeric indicating the significance of interest.
#' @param N Numeric indicating the length of data
#'
#' @return Numeric indicating the cutoff value
#' @noRd
.GeneralCutoff <- function(alpha, N){

  # Coeffs
  x <- -log(log(1-alpha)/-2)
  logN <- log(N)

  # Assumption 2.2
  a <- function(x){sqrt(2*log(x))}
  b <- function(x){2*log(x) + 1/2 * log(log(x)) - 1/2 * log(pi) }

  cutoff <- (x+b(logN)) / a(logN)

  cutoff
}


#' Vostrikova Cutoff
#'
#' Cufoff given in Theorem 2.1 with appropriate finite data approx dicussed in
#'     Vostrikova. See use in .computeCutoff.
#'
#' @param alpha Numeric indicating the significance of interest.
#' @param N Numeric indicating the length of data
#' @param nStart Numeric indicating starting value for trim.
#' @param nEnd Numeric indicating ending value for trim.
#'
#' @return Numeric indicating the cutoff value
#' @noRd
.VostrikovaCutoff <- function(alpha,N, nStart, nEnd){

  notDone <- TRUE
  low <- 0
  tmp <- 10
  while(is.na(tmp) || tmp>0){
    low <- low + 0.1
    tmp <- .findX(low, alpha = alpha, eAlpha=(N-(nEnd-nStart+1))/N)
  }
  ## P(sup U_st > x)
  cutoff <- stats::uniroot(f=.findX, lower=low,
                    upper=N, alpha=alpha, eAlpha=(N-(nEnd-nStart+1))/N)$root

  cutoff
}


#' Find X
#'
#' This (internal) function is used to compute x as needed in Vostrikova
#'     improvement of cutoff. See use in .VostrikovaCutoff.
#'
#' @param x Numeric x to consider
#' @param alpha Numeric indicating the significance of interest.
#' @param eAlpha Stardardized N^2 to handle values discussed 2.8
#'
#' @return Numeric indicating x value'
#' @noRd
.findX <- function(x, alpha, eAlpha){
  T <- 2*log((1-eAlpha)/eAlpha)
  alpha - ((x*exp(-x^2/2)/(2^(1/2)*sqrt(pi)))*(T-T/x^2+4/x^2))
}
