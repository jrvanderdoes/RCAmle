# Load data, format, and subset
UKcovid <- read.csv("./data-raw/UKcovid.csv")
UKcovid$date <- as.Date(UKcovid$date)
UKcovid <- UKcovid[UKcovid$date < as.Date('2022-1-1'), ]
UKcovid <- UKcovid[,c(2,4,5)]
colnames(UKcovid) <- c('nation','date','cases')

usethis::use_data(UKcovid, overwrite = TRUE)
