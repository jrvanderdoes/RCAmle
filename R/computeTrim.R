
#' Compute Trim
#'
#' This is just a function to control how much is trimmed
#'
#' @param requireDrops Integer indicating the number of require dropped. Ensure
#'     this is at least enough to estimate the parameters
#' @param type String of Start or End. This indicates if this is the start or
#'     ending trim
#' @param N (Optional) Integer of the length of data. Needed if type='End'
#'
#' @return Numeric indicating the index at which to start/stop examining in the
#'     data
#' @export
#'
#' @examples
#' computeTrim(10,'Start')
#' computeTrim(10,'End',100)
computeTrim <- function(requireDrops,type,N=NA){

  # if(is.na(mult)){
  if(type=='Start')
    return(requireDrops+1)
  else if(type=='End')
    return(N-requireDrops)
  # }
  # else{
  #   if(type=='Start'){
  #     if(!is.na(N))
  #       return(max(log(N),length(pars)*mult+1))
  #     return(length(pars)*mult+1)
  #   }else if(type=='End'){
  #     return(min(N-log(N), N-length(pars)*mult))
  #   }
  # }
}
