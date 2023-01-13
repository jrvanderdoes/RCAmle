## Load and organize the seperate files
folderLoc <- paste0("./data-raw/housing/")
cities <- c('Boston','Chicago','Denver','LasVegas','LosAngeles','Miami',
            'NewYork','SanDiego','SanFrancisco','WashingtonDC')

housing <- data.frame()#data.frame('city'=NA,'date'=NA,'index'=NA)

for(city in cities){
  cityData <- cbind(data.frame('city'=city),
                    as.data.frame(read.csv(paste0(folderLoc,city,'.txt'), sep='\t')))
  housing <- rbind(housing,cityData)
}

# Format result
colnames(housing) <- c('city','date','index')
housing$date <- as.Date(housing$date,'%d-%b-%Y')


usethis::use_data(housing, overwrite = TRUE)
