# Load Packages
library(RCAmle)
`%>%` <- magrittr::`%>%`
results <- list()


# The stationary region (Figure 3.1)
set.seed(12345)
data_sr <- simulateStationaryRegion(nSims = 10000000, silent=T)
data_sr[[3]]
results <- append(results, list('Fig3.1'=data_sr[[3]]))


# Empirical rejection frequencies under the null case (Table 3.1)
data_er <- simulationTable_NC(betas = c(0.5,0.75,1,1.05),
                              varProbRates = c(0.5, 0.5),
                              nSims=500, iterations = c(100, 200, 400, 800),
                              burnin = 1000, lowerEst=c(-Inf,0,10^-8,-Inf),
                              upperEst=c(Inf,Inf,Inf,Inf), alpha = 0.05,
                              errorTypes = c('Normal','Bernoulli','Exponential'),
                              CPLoc=0.5, seed=1234,trimAmt=10)

data_er$errType <- ifelse(data_er$errType=='Normal','N',
                          ifelse(data_er$errType=='Bernoulli','B','E'))
tmp <- data_er %>%
  dplyr::select(-c(Lower,Upper)) %>%
  tidyr::pivot_wider(names_from = c(iters,errType), values_from = Estim)

kableExtra::kable(tmp[,!(colnames(tmp)%in%c('beta'))], booktabs = TRUE) %>%
  kableExtra::pack_rows(index = table(tmp$beta))
results <- append(results, list('Tab3.1'=
                                  kableExtra::kable(tmp[,!(colnames(tmp)%in%c('beta'))], booktabs = TRUE) %>%
                                  kableExtra::pack_rows(index = table(tmp$beta))))


# The power curves (Figure 3.2)
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

# Figure 3.2, (upper left, upper right, lower left, lower right)
pc1[[1]]
pc2[[1]]
pc3[[1]]
pc4[[1]]
results <- append(results,
                  list('Fig3.2'=list(pc1[[1]],pc2[[1]],pc3[[1]],pc4[[1]])))


# The US housing figures (Figure 4.1)
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
                                       alpha=0.05,
                                       trimAmt = 10, silent = T )

  data_plot <- renderStackedPlot(trueData = data,
                     parCPData = result_Vost[,c(1:2,7:10)],
                     title = NULL,
                     subtitle = NULL,
                     varPlots = FALSE)

  ## Save data
  housingDataList <- append(housingDataList, list(list(result_Vost,data_plot)))
}

# Figure 4.1, (right, left)
housingDataList[[1]][[2]]
housingDataList[[2]][[2]]
results <- append(results,
                  list('Fig4.1'=list(housingDataList[[1]][[2]],
                                     housingDataList[[2]][[2]])))


# The UK Covid-19 figures (Figure 4.2)
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

# Figure 4.2 (upper left, upper right, lower left, lower right)
UKCovidList[[1]][[2]]
UKCovidList[[2]][[2]]
UKCovidList[[3]][[2]]
UKCovidList[[4]][[2]]
results <- append(results,
                  list('Fig4.2'=list(UKCovidList[[1]][[2]],
                                     UKCovidList[[2]][[2]],
                                     UKCovidList[[3]][[2]],
                                     UKCovidList[[4]][[2]])))


# Save RDS
saveRDS(results, 'results.rds')
