#' Simulate Stationary Region
#'
#' This simulates and draws the RCA(1) stationary region
#'
#' @param nSims (Optional) Integer indicating the number of simulations to use
#'     at each point
#' @param addEdge (Optional) Boolean indicating if the edge should be drawn
#' @param generationTrick (Optional) String that can indicate tricks to improve
#'     running time. Default is NA.
#'     \itemize{
#'         \item NA: No tricks, just run. This will take the most time
#'         \item theoretical: Rather than generating normals, use quantiles to
#'             get theoretical values
#'         \item reuse: Rather than generating unique normals at each point, use
#'             one set for all
#'         \item memory: Uses complete matrix multiplication rather than loops.
#'             Very fast, but likely will have memory issues on most machines
#'             for large nSims
#'     }
#' @param silent (Optional) Boolean indicating if output should be supressed (T).
#'     Default is FALSE.
#'
#' @return List with data.frame of date, (opt) data.frame of edge information,
#'     and ggplot2 plot.
#' @export
#'
#' @examples
#' dat <- simulateStationaryRegion()
simulateStationaryRegion <- function(nSims=10000, addEdge=TRUE,
                                     generationTrick=NA, silent=F){

  if(!silent) cat('-- Generate Region --\n')
  if(is.na(generationTrick)){
    stationaryDF <- .generateStationaryRegion_gen(
      seq(-1.5,1.5,0.01),
      seq(0, 4, 0.01),
      nSims,
      silent)
  }else if(generationTrick=='theoretical'){
    stationaryDF <- .generateStationaryRegion_theo(
      seq(-1.5,1.5,0.01),
      seq(0, 4, 0.01),
      nSims,
      silent)
  }else if(generationTrick=='memory'){
    stationaryDF <- .generateStationaryRegion_memory(
      seq(-1.5,1.5,0.01),
      seq(0, 4, 0.01),
      nSims)
  }else if(generationTrick=='reuse'){
    stationaryDF <- .generateStationaryRegion_reuse(
      seq(-1.5,1.5,0.01),
      seq(0, 4, 0.01),
      nSims,
      silent)
  } else{
    stop('Error: Please check generationTrick input in documentation and retry.')
  }

  ## Plots
  if(addEdge){
    if(!silent) cat('\n--- Compute Edge ---\n')
    edge <- .generateEdge(stationaryDF,silent)

    plot <- ggplot2::ggplot(stationaryDF, ggplot2::aes(x=var1,y=beta)) +
      ggplot2::geom_polygon(color='black', fill='gray',
                   size=0.75, data=edge) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size=22),
            plot.subtitle = ggplot2::element_text(hjust=0.5),
            legend.position="none",
            axis.title = ggplot2::element_text(size=18),
            axis.text = ggplot2::element_text(size=14)) +
      ggplot2::labs(y=bquote(beta),
           x=bquote('Var('~epsilon[1]~')'))

    return(list(stationaryDF, edge, plot))
  } else{

    plot <- ggplot2::ggplot(stationaryDF, ggplot2::aes(x=var1,y=beta)) +
      ggplot2::geom_tile(ggplot2::aes(fill=stable)) +
      ggplot2::theme_bw() +
      ggplot2::ggtitle("RCA Stationary",
              subtitle=paste0("Blue Stable, Red Unstable\n",
                              "Mean of ", nSims," trials"))+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
            plot.subtitle = ggplot2::element_text(hjust=0.5),
            legend.position="none")

    return(list(stationaryDF, plot))
  }
}


#' Generate Stationary Region - General
#'
#' This (internal) function computes stationarity of region using Monte Carlo
#'     simulation.
#'
#' @param bSeq Vector indicating sequence of beta values to be examined
#' @param epsSeq Vector indicating sequence of epsilon (1) values to be examined
#' @param nSims Integer indicating the number of simulation at each point
#' @param silent (Optional) Boolean indicating if output should be supressed (T).
#'     Default is FALSE.
#'
#' @return Data.frame indicating stationarity at each beta-epsilon1 value. The
#'     columns are beta (num), var1 (num), val (num), stable (bool).
#'
#' @examples
#' # This is an internal function. See use in simulateStationaryRegion.
.generateStationaryRegion_gen <- function(bSeq, epsSeq, nSims,silent=FALSE){

  stationaryDF <- expand.grid(bSeq, epsSeq)
  colnames(stationaryDF) <- c('beta','var1')
  stationaryDF$val <- NA

  N <- length(stationaryDF$beta)
  percent <- floor(N/10)
  if(!silent) st <- Sys.time()

  for(i in 1:N){
    ## Check for stationary
    #   stationary = E log|beta+e1|<0
    #     E log |pars[1]+N(0,sqrt(pars[2]))|
    stationaryDF$val[i] <-
      mean(
        log(
          abs(stationaryDF$beta[i] +
                rnorm(nSims,mean=0,
                      sd=sqrt(stationaryDF$var1[i])))
        )
      )

    ## Progress and Timing
    if(!(i %% percent) && !silent){
      indTime <- difftime(Sys.time(),st,units='mins')/i

      cat(paste0(i,'/',N,
                 ' (',round(indTime*(N-i),4),' min(s) remain)\n'))

    }
  }
  stationaryDF$stable <- stationaryDF$val<0

  stationaryDF
}


#' Generate Stationary Region - Memory
#'
#' This (internal) function computes stationarity of region using Monte Carlo
#'     simulation. This strictly uses matrices rather than loops. May have
#'     memory issues.
#'
#' @param bSeq Vector indicating sequence of beta values to be examined
#' @param epsSeq Vector indicating sequence of epsilon (1) values to be examined
#' @param nSims Integer indicating the number of simulation at each point
#'
#' @return Data.frame indicating stationarity at each beta-epsilon1 value. The
#'     columns are beta (num), var1 (num), val (num), stable (bool).
#'
#' @examples
#' # This is an internal function. See use in simulateStationaryRegion.
.generateStationaryRegion_memory <- function(bSeq, epsSeq, nSims){

  stationaryDF <- expand.grid(bSeq, epsSeq)
  colnames(stationaryDF) <- c('beta','var1')
  stationaryDF$val <- NA

  stationaryDF$val <-
    rowMeans( log(abs(stationaryDF$beta +
                        t(apply(t(sqrt(stationaryDF$var1)),
                                2,rnorm,n=nSims,mean=0))
    )))

  stationaryDF$stable <- stationaryDF$val<0

  stationaryDF
}


#' Generate Stationary Region - Reuse
#'
#' This (internal) function computes stationarity of region using Monte Carlo
#'     simulation. This reuses the same normals at each point to increase speed.
#'
#' @param bSeq Vector indicating sequence of beta values to be examined
#' @param epsSeq Vector indicating sequence of epsilon (1) values to be examined
#' @param nSims Integer indicating the number of simulation at each point
#' @param silent (Optional) Boolean indicating if output should be supressed (T).
#'     Default is FALSE.
#'
#' @return Data.frame indicating stationarity at each beta-epsilon1 value. The
#'     columns are beta (num), var1 (num), val (num), stable (bool).
#'
#' @examples
#' # This is an internal function. See use in simulateStationaryRegion.
.generateStationaryRegion_reuse <- function(bSeq, epsSeq, nSims, silent=FALSE){

  stationaryDF <- expand.grid(bSeq, epsSeq)
  colnames(stationaryDF) <- c('beta','var1')
  stationaryDF$val <- NA

  N <- length(stationaryDF$beta)
  percent <- floor(N/10)
  if(!silent) st <- Sys.time()

  norms <- rnorm(nSims)

  for(i in 1:N){
    ## Check for stationary
    #   stationary = E log|beta+e1|<0
    #     E log |pars[1]+N(0,sqrt(pars[2]))|
    stationaryDF$val[i] <-
      mean(
        log(
          abs(stationaryDF$beta[i] +
                sqrt(stationaryDF$var1[i]) * norms)
        )
      )

    ## Progress and Timing
    if(!(i %% percent) && !silent){
      indTime <- difftime(Sys.time(),st,units='mins')/i

      cat(paste0(i,'/',N,
                 ' (',round(indTime*(N-i),4),' min(s) remain)\n'))

    }
  }
  stationaryDF$stable <- stationaryDF$val<0

  stationaryDF
}


#' Generate Stationary Region - Theoretical
#'
#' This (internal) function computes stationarity of region using Monte Carlo
#'     simulation. This using quantiles rather than random computing MC
#'     simulations.
#'
#' @param bSeq Vector indicating sequence of beta values to be examined
#' @param epsSeq Vector indicating sequence of epsilon (1) values to be examined
#' @param nSims Integer indicating the number of simulation at each point
#' @param silent (Optional) Boolean indicating if output should be supressed (T).
#'     Default is FALSE.
#'
#' @return Data.frame indicating stationarity at each beta-epsilon1 value. The
#'     columns are beta (num), var1 (num), val (num), stable (bool).
#'
#' @examples
#' # This is an internal function. See use in simulateStationaryRegion.
.generateStationaryRegion_theo <- function(bSeq, epsSeq, nSims, silent=F){

  stationaryDF <- expand.grid(bSeq, epsSeq)
  colnames(stationaryDF) <- c('beta','var1')
  stationaryDF$val <- NA

  N <- length(stationaryDF$beta)
  percent <- floor(N/10)
  if(!silent) st <- Sys.time()

  norms <- qnorm(seq(1/(nSims+1),1-1/(nSims+1),1/(nSims+1)))

  for(i in 1:N){
    ## Check for stationary
    #   stationary = E log|beta+e1|<0
    #     E log |pars[1]+N(0,sqrt(pars[2]))|
    stationaryDF$val[i] <-
      mean(
        log(
          abs(stationaryDF$beta[i] +
                sqrt(stationaryDF$var1[i]) * norms)
        )
      )

    ## Progress and Timing
    if(!(i %% percent) && !silent){
      indTime <- difftime(Sys.time(),st,units='mins')/i

      cat(paste0(i,'/',N,
                 ' (',round(indTime*(N-i),4),' min(s) remain)\n'))

    }
  }
  stationaryDF$stable <- stationaryDF$val<0

  stationaryDF
}


#' Generate Edge
#'
#' This (internal) function finds the edges of the stationary region.
#'
#' @param data Data.frame indicating stationarity at each beta-epsilon1 value.
#'     The columns are beta (num), var1 (num), val (num), stable (bool).
#' @param silent (Optional) Boolean indicating if output should be supressed (T).
#'     Default is FALSE.
#'
#' @return Data,frame of data with column boundary (Bool) being bound indicating
#'     if the beta-var1 is a boundary position
#'
#' @examples
#' # This is an internal function. See use in simulateStationaryRegion.
.generateEdge <- function(data, silent=FALSE){
  data$boundary <- FALSE
  N <- length(data$beta)
  percent <- floor(N/10)

  if(!silent) st <- Sys.time()
  for(i in 1:N){

    if(data$stable[i]==1){

      currBeta <- data$beta[i]
      currVar1 <- data$var1[i]

      subsetData_upper <- data[data$beta > currBeta &
                                 data$var1==currVar1 &
                                 data$stable,]
      subsetData_lower <- data[data$beta < currBeta &
                                 data$var1==currVar1 &
                                 data$stable,]

      if(length(subsetData_upper$beta)==0 ||
         length(subsetData_lower$beta)==0)
        data$boundary[i] <- TRUE
    }

    ## Progress and Timing
    if(!(i %% percent) && !silent){
      indTime <- difftime(Sys.time(),st,units='mins')/i

      cat(paste0(i,'/',N,
                 ' (',round(indTime*(N-i),4),' min(s) remain)\n'))
    }
  }

  tmp <- data[data$boundary & data$beta>=0,]
  lineSort <- tmp[order(tmp$var1),]
  tmp <- data[data$boundary & data$beta<0,]
  lineSort <- rbind(lineSort,tmp[order(-tmp$var1),])

  return(lineSort)
}
