##Script to generate CDD models for sites 

library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(suncalc)
library(rnoaa)
library(doParallel)

n.cores <- 20
registerDoParallel(cores=n.cores)
source('phenologyModel.R')

siteData <- read.csv("data/phenologyForecastSites.csv",header=TRUE)
endDate <- as.Date("2019-12-31")
dataDirectory <- "data/"
baseTemp <- 20
vars <- "noCov"

#foreach(i=1:nrow(siteData)) %dopar% {
for(i in 1:nrow(siteData)){
  lat <- as.numeric(siteData[i,2])
  long <- as.numeric(siteData[i,3])
  startDate <- (as.Date(siteData[i,7]))
  URL <- as.character(siteData$URL[i])
  URL2 <- as.character(siteData$URL2[i])
  URL3 <- as.character(siteData$URL3[i])
  if(nchar(URL2)>0){
    URL <- c(URL,URL2)
    if(nchar(URL3)>0){
      URL <- c(URL,URL3)
    }
  }
  TZ <- as.numeric(siteData[i,6])
  URLs <- URL
  load(file=paste("PhenologyForecastData/",siteName,"_phenopixOutputs.RData",sep=""))
  fittedDat=allDat
  ERA5dataFolder <- paste("/projectnb/dietzelab/kiwheel/ERA5/Data/",siteName,"/",sep="")
  phenologyCalibration_Autumn(siteName=siteName,URLs=URLs,lat=lat,long=long,dataDirectory="PhenologyForecastData/",
                              startDate=startDate,endDate=endDate,baseTemp=baseTemp,
                              ERA5dataFolder = ERA5dataFolder,vars=vars,TZ=TZ,calDIC=FALSE,fittedDat=allDat,
                              splitYears="first")
}


