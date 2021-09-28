##Script to create consolidated data objects (dataFinal)

library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(suncalc)
library(rnoaa)
library(doParallel)
library(ncdf4)
source('load_ERA5.R')

n.cores <- 8

#register the cores.
registerDoParallel(cores=n.cores)

dataDirectory <- "data/"
siteData <- read.csv('/projectnb/dietzelab/kiwheel/chlorophyllCycling/allPhenocamDBsitesComplete.csv',header=TRUE)

#for(s in 1:nrow(siteData)){
foreach(s=1:nrow(siteData)) %dopar% {
  siteName <- as.character(siteData$siteName[s])
  lat <- as.numeric(siteData[s,2])
  long <- as.numeric(siteData[s,3])
  startDate <- (as.Date(siteData[s,7]))
  endDate <- as.character(siteData$endDate[s])
  URL <- as.character(siteData$URL[s])
  URL2 <- as.character(siteData$URL2[s])
  URL3 <- as.character(siteData$URL3[s])
  if(!is.na(URL2)){
    URL <- c(URL,URL2)
    if(!is.na(URL3)){
      URL <- c(URL,URL3)
    }
  }
  TZ <- as.numeric(siteData[s,6])
  URLs <- URL
  load(file=paste(dataDirectory,siteName,"_phenopixOutputs.RData",sep=""))
  fittedDat=allDat
  ERA5dataFolder <- paste("/projectnb/dietzelab/kiwheel/ERA5/Data/",siteName,"/",sep="")
  
  phenoData <- matrix(nrow=0,ncol=32)
  print(URLs[1])
  for(u in 1:length(URLs)){
    phenoDataSub <- download.phenocam(URLs[u])
    phenoData <- rbind(phenoData,phenoDataSub)
  }
  
  ##Order and remove duplicate PC data
  phenoData2 <- phenoData[order(phenoData$date),]
  phenoData3 <- phenoData2[!duplicated(phenoData2$date),]
  phenoData <- phenoData3
  
  phenoData <- phenoData[phenoData$date<endDate,]
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
  calFileName <- paste0(siteName,"_",startDate,"_",endDate,"_era5TemperatureMembers.nc")
  datTairEns <- load_ERA5(ERA5dataFolder=ERA5dataFolder,calFileName=calFileName,TZ_offset=TZ,variable="Tair")
  
  TairMu <- apply(X=datTairEns,MARGIN=2,FUN=mean)
  TairPrec <- 1/apply(X=datTairEns,MARGIN=2,FUN=var)
  dat2$TairMu <- TairMu 
  dat2$TairPrec<- TairPrec
  
  dayLengths <- numeric()
  
  for(d in 1:length(days)){
    suntimes <- getSunlightTimes(date=days[d],
                                 lat=lat,lon=long,keep=c("nauticalDawn","nauticalDusk"),
                                 tz = "GMT") #GMT because I only care about difference
    dayLengths <- c(dayLengths,as.numeric(suntimes$nauticalDusk-suntimes$nauticalDawn))
  }
  
  dat2$D <- dayLengths
  ICsdat <- dat2[as.numeric(format(dat2$dates,"%j"))%in% seq(203,212),]
  dat2 <- dat2[as.numeric(format(dat2$dates,"%j"))%in% seq(213,365),]
  
  nrowNum <- 365-212
  p <- matrix(nrow=nrowNum,ncol=0)
  TairMu <- matrix(nrow=nrowNum,ncol=0)
  D <- matrix(nrow=nrowNum,ncol=0)
  ICs <- matrix(nrow=10,ncol=0)
  TairPrec <- matrix(nrow=nrowNum,ncol=0)
  valNum <- 0
  days2 <- matrix(nrow=nrowNum,ncol=0)
  
  finalYrs <- numeric()
  sofs <- numeric()
  for(i in (lubridate::year(as.Date(dat2$dates[1]))):lubridate::year(as.Date(dat2$dates[length(dat2$dates)]))){
    subDat <- dat2[lubridate::year(as.Date(dat2$dates))==i,]
    valNum <- valNum + 1
    Low <- fittedDat[valNum,'Low']
    High <- fittedDat[valNum,'High']
    if(!is.na(Low)){
      newCol <- scales::rescale(subDat$p,to=c(0,1),from=c(Low,High))
      p <- cbind(p,newCol)
      ICs <- cbind(ICs,scales::rescale(ICsdat[lubridate::year(as.Date(ICsdat$dates))==i,]$p,from=c(Low,High)))
      days2 <- cbind(days2,as.Date(subDat$dates))
      finalYrs <- c(finalYrs,i)
      sofs <- c(sofs,(fittedDat[valNum,'FallStartDay']-212)) ######Change for start if needed
      TairMu <- cbind(TairMu,subDat$TairMu)
      D <- cbind(D,subDat$D)
      TairPrec <- cbind(TairPrec,subDat$TairPrec)
    }
  }
  p[p<0] <- 0
  p[p>0.999] <- 0.999
  ICs[ICs<0] <- 0
  ICs[ICs>0.999] <- 0.999
  
  dataFinal <- list(p=p,years=finalYrs,sofMean=mean(sofs))
  dataFinal$n <- nrowNum
  dataFinal$N <- ncol(dataFinal$p)
  
  
  x1a <- numeric()
  x1b <- numeric()
  for(yr in 1:dataFinal$N){
    mu <- mean(ICs[,yr],na.rm=TRUE)
    vr <- var(ICs[,yr],na.rm = TRUE)
    x1a <- c(x1a,(mu**2-mu**3-mu*vr)/(vr))
    x1b <- c(x1b,(mu-2*mu**2+mu**3-vr+mu*vr)/(vr))
  }
  
  dataFinal$x1.a <- x1a
  dataFinal$x1.b <- x1b
  
  dataFinal$TairMu <- TairMu
  dataFinal$TairPrec <- TairPrec
  dataFinal$D <- D
  
  save(dataFinal,file=paste0(dataDirectory,siteName,"_dataFinal.RData"))
  
}