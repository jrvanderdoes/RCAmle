
#' Likelihood
#'
#' This (internal) function computes the liklihood, defined as 2.2 in the paper
#'
#' @param y Vector of numerics for the data
#' @param u Vector of numerics for parameters, c(beta, sigma1^2, sigma2^2)
#'
#' @return Vector of numerics indicating the liklihoods at each point (removing
#'     uncomputable ones)
#' @noRd
.l <- function(y,u){

  l <- rep(NA,length(y))

  for(i in 2:length(y)){
    l[i] <- -1/2 * ( log(u[2] * y[i-1]^2 + u[3]) +
                       (y[i] - u[1] * y[i-1])^2 / (u[2]*y[i-1]^2 + u[3]) )
  }

  l[-1]
}


#' Log Liklihood
#'
#' This (internal) function computes the log liklihood, defined as 2.1 in the paper
#'
#' @param u Vector of numerics for parameters, c(beta1, sigma1^2, sigma2^2, beta2)
#' @param k Integer indicating the value to stop using beta1 (i.e. the change point)
#' @param y Vector of numerics for the data
#' @param N (Optional) Integer indicating the length of y. Value added to
#'     increase speed rather than calling length(y) each call.
#'
#' @return Numeric of log liklihood with change point k
#' @noRd
.L <- function(u, k, y, N=length(y)){

  ## Second goes from k since drops first (so k+1:N)
  L <- sum(.l(y[1:k], c(u[1], u[2], u[3]))) +
    sum(.l(y[k:N], c(u[4],u[2],u[3])))

  L
}


#' (Array) Log Liklihood
#'
#' This (internal) function computes the log liklihood, defined as 2.1 in the
#'     paper. This is built to better fit for array calls.
#'
#' @param u1 Numeric for parameter 1, beta1
#' @param u2 Numeric for parameter 2, sigma1^2
#' @param u3 Numeric for parameter 3, sigma2^2
#' @param u4 Numeric for parameter 4, beta2
#' @param k Integer indicating the value to stop using beta1 (i.e. the change point)
#' @param y Vector of numerics for the data
#' @param N (Optional) Integer indicating the length of y. Value added to
#'     increase speed rather than calling length(y) each call.
#'
#' @return Numeric of log liklihood with change point k
#' @noRd
.L_array <- function(u1, u2, u3, u4, k, y, N = length(y)){

  y <- unlist(y)

  ## Second goes from k since drops first (so k+1:N)
  L <- sum(.l(y[1:k], c(u1, u2, u3))) +
    sum(.l(y[k:N], c(u4, u2, u3)))

  L
}
