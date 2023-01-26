---
title: "Accompanying Paper Code"
#output: rmarkdown::html_vignette
output: pdf_document
vignette: >
  %\VignetteIndexEntry{paperCode}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---





This vignette is a companion to the paper `The maximally selected likelihood 
ratio test in random coefficient models.' The code to build all figures in that
paper are presented. Discussions on these figures are conducted in the paper and
are not addressed in this data supplement.

The stationary region was built using the following code.

```r
set.seed(12345)
data_sr <- simulateStationaryRegion(nSims = 10000000, silent=T)
```

The results of that code includes the following figure.



```r
data_sr[[3]]
```

![Simulation stationary region using Gaussian error based on $E\log |\beta + \epsilon_{0,1} |$](paperCode_files/figure-latex/plotStationaryRegion-1.pdf) 


The empirical rejection frequencies under the null case was generated using the following code.

```r
data_er <- simulationTable_NC(betas = c(0.5,0.75,1,1.05), 
                              varProbRates = c(0.5, 0.5),
                              nSims=500, iterations = c(100, 200, 400, 800), 
                              burnin = 1000, lowerEst=c(-Inf,0,10^-8,-Inf), 
                              upperEst=c(Inf,Inf,Inf,Inf), alpha = 0.05,
                              errorTypes = c('Normal','Bernoulli','Exponential'),
                              CPLoc=0.5, seed=1234,trimAmt=10)
```

The results of that code is the following data.



```r
data_er$errType <- ifelse(data_er$errType=='Normal','N',
                          ifelse(data_er$errType=='Bernoulli','B','E'))
tmp <- data_er %>%
  dplyr::select(-c(Lower,Upper)) %>%
  pivot_wider(names_from = c(iters,errType), values_from = Estim)

kable(tmp[,!(colnames(tmp)%in%c('beta'))], booktabs = TRUE) %>%
  pack_rows(index = table(tmp$beta))
```


\begin{tabular}{lrrrrrrrrrrrr}
\toprule
type & 100\_N & 200\_N & 400\_N & 800\_N & 100\_B & 200\_B & 400\_B & 800\_B & 100\_E & 200\_E & 400\_E & 800\_E\\
\midrule
\addlinespace[0.3em]
\multicolumn{13}{l}{\textbf{0.5}}\\
\hspace{1em}MLE & 0.002 & 0.012 & 0.004 & 0.010 & 0.002 & 0.006 & 0.006 & 0.004 & 0.006 & 0.010 & 0.028 & 0.016\\
\hspace{1em}Vost & 0.026 & 0.046 & 0.024 & 0.032 & 0.032 & 0.032 & 0.034 & 0.028 & 0.028 & 0.026 & 0.042 & 0.052\\
\hspace{1em}WLS & 0.004 & 0.012 & 0.002 & 0.018 & 0.012 & 0.012 & 0.012 & 0.010 & 0.004 & 0.006 & 0.014 & 0.016\\
\addlinespace[0.3em]
\multicolumn{13}{l}{\textbf{0.75}}\\
\hspace{1em}MLE & 0.002 & 0.010 & 0.002 & 0.010 & 0.010 & 0.006 & 0.004 & 0.004 & 0.004 & 0.008 & 0.022 & 0.026\\
\hspace{1em}Vost & 0.030 & 0.028 & 0.028 & 0.036 & 0.044 & 0.024 & 0.032 & 0.026 & 0.026 & 0.024 & 0.048 & 0.052\\
\hspace{1em}WLS & 0.020 & 0.032 & 0.026 & 0.034 & 0.034 & 0.024 & 0.056 & 0.040 & 0.006 & 0.004 & 0.018 & 0.024\\
\addlinespace[0.3em]
\multicolumn{13}{l}{\textbf{1}}\\
\hspace{1em}MLE & 0.008 & 0.010 & 0.002 & 0.006 & 0.006 & 0.010 & 0.008 & 0.012 & 0.006 & 0.008 & 0.016 & 0.020\\
\hspace{1em}Vost & 0.022 & 0.040 & 0.032 & 0.036 & 0.048 & 0.032 & 0.036 & 0.038 & 0.030 & 0.036 & 0.054 & 0.040\\
\hspace{1em}WLS & 0.030 & 0.054 & 0.058 & 0.058 & 0.076 & 0.098 & 0.138 & 0.142 & 0.002 & 0.004 & 0.018 & 0.016\\
\addlinespace[0.3em]
\multicolumn{13}{l}{\textbf{1.05}}\\
\hspace{1em}MLE & 0.002 & 0.016 & 0.012 & 0.012 & 0.010 & 0.010 & 0.008 & 0.008 & 0.004 & 0.010 & 0.020 & 0.024\\
\hspace{1em}Vost & 0.024 & 0.052 & 0.032 & 0.044 & 0.046 & 0.038 & 0.038 & 0.036 & 0.028 & 0.036 & 0.058 & 0.048\\
\hspace{1em}WLS & 0.026 & 0.078 & 0.058 & 0.066 & 0.066 & 0.090 & 0.136 & 0.150 & 0.004 & 0.006 & 0.014 & 0.016\\
\bottomrule
\end{tabular}


The power curves were built using the following code.

```r
tmp <- c(1, 0.9,0.75, seq(0.6,0.3,-0.1),0.25,0.1)
U4seq <- c(tmp, 0, rev(-tmp))
set.seed(1234)
pc1 <- simulatePowerCurve(Us3 = c(0, 0.5, 0.5), U4s = U4seq,
                     nSims = 500, lowerEst = c(-Inf,0,10^-8,-Inf), 
                     upperEst = rep(Inf,4), alpha = 0.05,
                     burnin = 1000, k = 0.5 * 400, N = 400, 
                     errorType = 'Normal', trimAmt=10,
                     silent=T)

set.seed(1234)
pc2 <- simulatePowerCurve(Us3 = c(0, 0.5, 0.5), U4s = U4seq,
                     nSims = 500, lowerEst = c(-Inf,0,10^-8,-Inf), 
                     upperEst = rep(Inf,4), alpha = 0.05,
                     burnin = 1000, k = 0.9 * 400, N = 400, 
                     errorType = 'Normal', trimAmt=10,
                     silent=T)

U4seq <- c(1.05,seq(1,0.6,-0.1), 0.5, seq(0.4,0,-0.1),-0.05)
set.seed(1234)
pc3 <- simulatePowerCurve(Us3 = c(0.5, 0.5, 0.5), U4s = U4seq,
                     nSims = 500, lowerEst = c(-Inf,0,10^-8,-Inf), 
                     upperEst = rep(Inf,4), alpha = 0.05,
                     burnin = 1000, k = 0.5 * 400, N = 400, 
                     errorType = 'Normal', trimAmt=10,
                     silent=T)
set.seed(1234)
pc4 <- simulatePowerCurve(Us3 = c(0.5, 0.5, 0.5), U4s = U4seq,
                     nSims = 500, lowerEst = c(-Inf,0,10^-8,-Inf), 
                     upperEst = rep(Inf,4), alpha = 0.05,
                     burnin = 1000, k = 0.5 * 400, N = 400,
                     errorType = 'Bernoulli', trimAmt=10,
                     silent=T)
```

The resulting power figures are as follows (See the paper for details on each of the figures).



```r
pc1[[1]]
```

![Power curves as a function of the change in RCA coefficient.](paperCode_files/figure-latex/plotPowerCurves-1.pdf) 

```r
pc2[[1]]
```

![Power curves as a function of the change in RCA coefficient.](paperCode_files/figure-latex/plotPowerCurves-2.pdf) 

```r
pc3[[1]]
```

![Power curves as a function of the change in RCA coefficient.](paperCode_files/figure-latex/plotPowerCurves-3.pdf) 

```r
pc4[[1]]
```

![Power curves as a function of the change in RCA coefficient.](paperCode_files/figure-latex/plotPowerCurves-4.pdf) 


The US housing figures were built using the following code.

```r
housingDataList <- list()
for(city in c('Boston','LosAngeles')){
  
  ## Load data
  data <- housing[housing$city==city,c(2,3)]
  data <- data[order(data$date),]
  
  ## Detect CP
  set.seed(12345)
  result_Vost <- binarySegmentationCPDetection(fullData=data, 
                                       method='Vostrikova',
                                       lower=c(-Inf, 0, 10^-8, -Inf),
                                       upper=c(Inf, Inf, Inf, Inf),
                                       alpha=0.05, nStart=NA, nEnd=NA,
                                       trimAmt = 10, silent = T )
  
  data_plot <- renderStackedPlot(trueData = data, 
                     parCPData = result_Vost[,c(1:2,7:10)],
                     title = NULL,
                     subtitle = NULL,
                     varPlots = FALSE)
  
  ## Save data
  housingDataList <- append(housingDataList, list(list(result_Vost,data_plot)))
}
```

The results were as follows.



```r
housingDataList[[1]][[2]]
```

![Daily US housing price indices in Boston.](paperCode_files/figure-latex/plotHousingData_Boston-1.pdf) 

```r
housingDataList[[2]][[2]]
```

![Daily US housing price indices in Los Angeles.](paperCode_files/figure-latex/plotHousingData_LA-1.pdf) 


The UK covid-19 figures were built using the subsequent code.

```r
UKCovidList <- list()
for(nation in c('England','Northern Ireland','Scotland','Wales')){
  
  ## Load data
  data <- UKcovid[UKcovid$nation==nation,c(2,3)]
  data <- data[order(data$date),]
  
  ## Detect CP
  set.seed(12345)
  result_Vost <- binarySegmentationCPDetection(fullData=data, 
                                       method='Vostrikova',
                                       lower=c(-Inf, 0, 10^-8, -Inf),
                                       upper=c(Inf, Inf, Inf, Inf),
                                       alpha=0.05,
                                       trimAmt = 10, silent = T )
  data_plot <- renderStackedPlot(trueData = data, 
                     parCPData = result_Vost[,c(1:2,7:10)],
                     title = NULL,
                     subtitle = NULL,
                     varPlots = FALSE)
  
  ## Save data
  UKCovidList <- append(UKCovidList, list(list(result_Vost,data_plot)))
}
```

With plots produced as follows.



```r
UKCovidList[[1]][[2]]
```

![\label{fig:UKCovid}Daily Covid-19 patients hospitalised in England.](paperCode_files/figure-latex/plotUKCovidData_England-1.pdf) 

```r
UKCovidList[[2]][[2]]
```

![\label{fig:UKCovid}Daily Covid-19 patients hospitalised in Northern Ireland.](paperCode_files/figure-latex/plotUKCovidData_NIreland-1.pdf) 

```r
UKCovidList[[3]][[2]]
```

![\label{fig:UKCovid}Daily Covid-19 patients hospitalised in Scotland.](paperCode_files/figure-latex/plotUKCovidData_Scotland-1.pdf) 

```r
UKCovidList[[4]][[2]]
```

![\label{fig:UKCovid}Daily Covid-19 patients hospitalised in Wales.](paperCode_files/figure-latex/plotUKCovidData_Wales-1.pdf) 
