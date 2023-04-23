
#' Compute Amount of Data to Trim
#'
#' This is just a function to control how much is trimmed
#'
#' @param requireDrops Integer indicating the number of require dropped. Ensure
#'     this is at least enough to estimate the parameters. If NA, use log(N)
#'     for trim (min of 4)
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
  if(type=='Start'){
    if(is.na(requireDrops))
      return(max(4,ceiling(log(N)))+1)
    return(requireDrops+1)
  }else if(type=='End'){
    if(is.na(requireDrops))
      return(min(N-4,N-ceiling(log(N))))
    return(N-requireDrops)
  }
}
