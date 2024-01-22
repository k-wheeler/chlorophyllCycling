library(phenopix)
library(zoo)
library(PhenoForecast)


estimateDatesAndScale <- function(p,siteName,year){
  newP <- zoo(na.approx(p))
  #out <- KlostermanFit(newP)
  out <- ElmoreFit(newP)
  output <- c(year,out$fit$sf)
  x=out$fit$predicted
  dts <- PhenoDeriv(x=x,fit=out$fit,uncert=TRUE)
  output <- c(output,dts[4])
  xd <- c(NA,diff(x))
  xd2 <- c(NA,diff(xd))
  #xd3 <- c(NA,diff(xd2))
  sof <- 200+which.min(xd2[200:length(xd2)])-3
  #LastMax <- which.max(xd2[200:length(xd2)])
  #sof <- 200+which.min(xd2[200:(200+LastMax)])-3
  
  output <- c(output,sof)
  output <- c(output,
              which(out$fit$predicted[output[4]:365]<(output[3]-(output[3]-output[2])/2))[1])
  PhenoPlot(out,metrics="GCC",main=paste(siteName,year))
  points(newP,col="green",pch=20)
  abline(v=output[c(4,5)],col="red")
  return(output)
}
checkForPoorFits <- function(siteName,allDat,badSiteYears){ #Based on visual expection of years
  yrs <- na.omit(as.numeric(badSiteYears[badSiteYears[,1]==siteName,][2:ncol(badSiteYears)]))
  allDat[allDat[,'Year'] %in% yrs,2:6] <- c(NA,NA,NA,NA,NA)
  return(allDat)
}

calculateFits <- function(siteName,URLs,startDate,endDate,year){
  ###Download PhenoCam data and format
  phenoData <- matrix(nrow=0,ncol=32)
  print(URLs[1])
  for(u in 1:length(URLs)){
    print(URLs[u])
    phenoDataSub <- download.phenocam(URLs[u])
    phenoData <- rbind(phenoData,phenoDataSub)
  }
  
  ##Order and remove duplicate PC data
  phenoData2 <- phenoData[order(phenoData$date),]
  phenoData3 <- phenoData2[!duplicated(phenoData2$date),]
  phenoData <- phenoData3
  
  phenoData <- phenoData[phenoData$date<=endDate,]
  p.old <- phenoData$gcc_90
  time.old <-  as.Date(phenoData$date)
  
  days <- seq(as.Date(startDate),(as.Date(endDate)),"day")
  p <- rep(NA,length(days))
  
  for(i in 1:length(p.old)){
    p[which(days==time.old[i])] <- p.old[i]
  }
  
  months <- lubridate::month(days)
  years <- lubridate::year(days)
  
  dat2 <- data.frame(dates=days,years=years,months=months,p=p)
  
  dat2 <- dat2[as.numeric(format(dat2$dates,"%j"))%in% seq(1,365),] #Remove final day in leap years
  allDat <- matrix(nrow=0,ncol=6)
  
  for(i in (lubridate::year(as.Date(dat2$dates[1]))):lubridate::year(as.Date(dat2$dates[length(dat2$dates)]))){##I know this includes the forecasted stuff, but it shouldn't really matter because of the JAGS model setup
    print(i)
    subDat <- dat2[lubridate::year(as.Date(dat2$dates))==i,]
    if(length(na.omit(subDat$p))>1){
      allDat <- rbind(allDat,estimateDatesAndScale(p=subDat$p,siteName,year=i))
    }
  }
  colnames(allDat) <- c("Year","Low","High","PeakDay","FallStartDay","perc50Day")
  return(allDat)
}

#siteData <- read.csv('/projectnb/dietzelab/kiwheel/chlorophyllCycling/allPhenocamDBsitesComplete.csv',header=TRUE)
siteData <- read.csv('allPhenocamDBsitesComplete.csv',header=TRUE)
badSiteYears <- read.csv('phenocamSitesBadYears.csv',header = FALSE)
dataDirectory <- "data/"
#endDate <- as.Date("2020-12-31")

pdf(file="ChlorophyllCycling_AutumnForecastPreviousCurveElmoreFits.pdf",height=5,width=8)
for(i in 1:nrow(siteData)){
#for(i in 1:10){
  #for(i in 1:5){
  siteName <- as.character(siteData[i,1])
  print(siteName)
  
  startDate <- (as.Date(siteData[i,7]))
  endDate <- as.Date(siteData$endDate[i])
  URL <- as.character(siteData$URL[i])
  URL2 <- as.character(siteData$URL2[i])
  URL3 <- as.character(siteData$URL3[i])
  if(!is.na(URL2)){
    URL <- c(URL,URL2)
    if(!is.na(URL3)){
      URL <- c(URL,URL3)
    }
  }
  URLs <- URL
  allDat <- calculateFits(siteName,URLs,startDate,endDate)
  allDat <- checkForPoorFits(allDat=allDat,siteName=siteName,badSiteYears=badSiteYears)
  save(allDat,file=paste(dataDirectory,siteName,"_phenopixOutputs2.RData",sep=""))
}
dev.off()
