##Script to Populate the PhenoCam data table:
library(PhenologyBayesModeling)
valid_url <- function(url_in,t=2){
  con <- url(url_in)
  check <- suppressWarnings(try(open.connection(con,open="rt",timeout=t),silent=T)[1])
  suppressWarnings(try(close.connection(con),silent=T))
  ifelse(is.null(check),TRUE,FALSE)
}

dat <- read.csv("allPhenocamDBsites.csv",header=TRUE)
dat$startDate <- as.Date(1:nrow(dat), origin=Sys.Date())
dat$endDate <- as.Date(1:nrow(dat), origin=Sys.Date())

for(i in 1:nrow(dat)){
#for(i in 1:15){
  siteName <- as.character(dat$siteName[i])
  print(siteName)
  
  URL_table <- data.frame(matrix(nrow=8,ncol=4))
  colnames(URL_table) <- c("URL","exists","startDate","endDate")
  URL_table$URL <- c(paste0("http://phenocam.sr.unh.edu/data/archive/",siteName,"/ROI/",
                            siteName,"_DB_0001_1day.csv"),
                     paste0("http://phenocam.sr.unh.edu/data/archive/",siteName,"/ROI/",
                            siteName,"_DB_0002_1day.csv"),
                     paste0("http://phenocam.sr.unh.edu/data/archive/",siteName,"/ROI/",
                            siteName,"_DB_0003_1day.csv"),
                     paste0("http://phenocam.sr.unh.edu/data/archive/",siteName,"/ROI/",
                            siteName,"_DB_0004_1day.csv"),
                     paste0("http://phenocam.sr.unh.edu/data/archive/",siteName,"/ROI/",
                            siteName,"_DB_1000_1day.csv"),
                     paste0("http://phenocam.sr.unh.edu/data/archive/",siteName,"/ROI/",
                            siteName,"_DB_2000_1day.csv"),
                     paste0("http://phenocam.sr.unh.edu/data/archive/",siteName,"/ROI/",
                            siteName,"_DB_3000_1day.csv"),
                     paste0("http://phenocam.sr.unh.edu/data/archive/",siteName,"/ROI/",
                            siteName,"_DB_4000_1day.csv"))
  
  URL_table$startDate <- as.Date(1:8, origin=Sys.Date())
  URL_table$endDate <- as.Date(1:8, origin=Sys.Date())
  for(u in 1:nrow(URL_table)){
    URL_table$exists[u] <- valid_url(url_in=URL_table$URL[u])
    URL_table$exists[u] <- URL_table$exists[u] || valid_url(url_in=URL_table$URL[u]) ##Tries three times to make sure 
    URL_table$exists[u] <- URL_table$exists[u] || valid_url(url_in=URL_table$URL[u])
    if(URL_table$exists[u]){
      if(is.na(dat$TZ[i])){
        phenoDat <- read.csv(URL_table$URL[u])
        dat$TZ[i] <- strsplit(as.character(phenoDat[9,1])," ")[[1]][4]
      }
      phenoDat <- read.csv(URL_table$URL[u],skip = 22)
      URL_table$startDate[u] <- as.Date(range(as.Date(phenoDat$date)))[1]
      URL_table$endDate[u] <- as.Date(range(as.Date(phenoDat$date)))[2]
    }else{
      URL_table$startDate[u] <- NA
      URL_table$endDate[u] <- NA
    }
  }
  if(sum(URL_table$exists)>0){
  maxLength <- which.max(URL_table$endDate - URL_table$startDate)
  minVals <- which(URL_table$startDate==min(URL_table$startDate,na.rm = TRUE))
  maxVals <- which(URL_table$endDate==max(URL_table$endDate,na.rm = TRUE))
  if(maxLength %in% minVals){
    finalURLs <- URL_table$URL[maxLength]
  }else{
    finalURLs <- c(URL_table$URL[minVals[1]],
                   URL_table$URL[maxLength])
  }
  if(!(maxLength %in% maxVals)){
    finalURLs <- c(finalURLs,URL_table$URL[maxVals[1]])
  }
  startDate <- min(URL_table$startDate,na.rm=TRUE)
  if(!(lubridate::month(startDate)==1 & lubridate::day(startDate)==1)){
    startYear <- lubridate::year(startDate) + 1
    startDate <- as.Date(paste0(startYear,"-01-01"))
  }
  dat$startDate[i] <- startDate
  endDate <- max(URL_table$endDate,na.rm = TRUE) 
  if(!(lubridate::month(endDate)==12 & lubridate::day(endDate)==31)){
    endYear <- lubridate::year(endDate) - 1
    endDate <- as.Date(paste0(endYear,"-12-31"))
  }
  dat$endDate[i] <- endDate
  dat$URL[i] <- finalURLs[1]
  if(length(finalURLs>1)){
    dat$URL2[i] <- finalURLs[2]
    if(length(finalURLs>2)){
      dat$URL3[i] <- finalURLs[3]
    }
  }

  }else{
    dat$startDate[i] <- NA
    dat$endDate[i] <- NA
  }

  
}

dat2 <- dat
dat <- dat[!is.na(dat$URL),]
for(i in 1:nrow(dat)){
  if(dat$startDate[i]>=dat$endDate[i]){
    dat$startDate[i] <- NA
    dat$endDate[i] <- NA
    dat$URL[i] <- NA
  }
}
dat <- dat[!is.na(dat$URL),]

dat[dat$siteName=="harvardlph",]$startDate <- as.Date("2016-01-01")
dat[dat$siteName=="howland2",]$startDate <- as.Date("2011-01-01")

write.csv(dat,file="allPhenocamDBsitesComplete.csv",row.names = FALSE,quote=FALSE)

          
