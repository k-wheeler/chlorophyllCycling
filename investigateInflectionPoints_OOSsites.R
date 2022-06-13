library('rjags')
library('runjags')
library('ecoforecastR')
library('RColorBrewer')
library(doParallel)
#Create forecast step model:
#' Basic logistic forecast step
#'
#' @param IC Initial conditions
#' @param b0 The parameter b0
#' @param b1 The parameter b1
#' @param b2 The parameter b2
#' @param b3 The parameter b3
#' @param b4 The parameter b4
#' @param b5 The parameter b5
#' @param Q Process error (default = 0 for deterministic runs)
#' @param n Size of Monte Carlo ensemble
#' @param NT number of days for forecast
#' @param Tair
#' @param D
#' @import rjags
#' @import runjags
#' @import ecoforecastR
#'
#' @return
#' @export
#'
#' @examples
forecastStep <- function(IC,b0,b1,b2,b3,b4,Q=0,n,NT,Tair,D){
  x <- matrix(NA,n,NT)
  if(length(IC)==1){
    Xprev <- rep(IC,n)
  }else{
    Xprev <- IC
  }
  
  for(t in 1:NT){
    bd <- Xprev + b4 * Xprev #b0 + (b1 * Xprev) + (b2 * Xprev **2) 
    syn <- b0 + b1 * Tair[t] + b2 * D[t] + b3 * Tair[t] * D[t] 
    if(length(syn)==1){
      syn <- rep(syn,n)
    }
    
    xNew <- numeric()
    mu <- rep(NA,n)
    for(i in 1:n){
      mu[i] <- bd[i] + max(syn[i],0)
      mu[i] <- max(0,min(mu[i],IC[min(i,length(IC))]))
      
      if(length(Q)>1){
        xNew <- c(xNew,rnorm(1,mu[i],Q[i]))
      }else{
        xNew <- c(xNew,rnorm(1,mu[i],Q))
      }
    }
    x[,t] <- xNew
    Xprev <- x[,t]
  }
  return(x)
}

calculateStart3 <- function(ys,avgNum=10){
  avgDiffs <- rep(0,avgNum)
  for(t in (avgNum + 1):(length(ys)-avgNum)){
    prevAvg <- mean(diff(ys[1:(t-1)]),na.rm = TRUE)
    newAvg <- mean(diff(ys[t:(t+avgNum)]),na.rm = TRUE)
    avgDiffs <- c(avgDiffs,(newAvg-prevAvg))
  }
  avgDiffs <- c(avgDiffs,rep(0,avgNum))
  return(avgDiffs)
}

dataDirectory <- "data/"
siteData <- read.csv('/projectnb/dietzelab/kiwheel/chlorophyllCycling/allPhenocamDBsitesComplete.csv',header=TRUE)
#siteName <- "lacclair"
#calSite <- "lacclair"
n.cores <- 24
registerDoParallel(cores=n.cores)

calSites <- c("harvard","umichbiological","bostoncommon","coweeta","howland2",
           "morganmonroe","missouriozarks","queens","dukehw","lacclair","bbc1","NEON.D08.DELA.DP1.00033",
           "bartlettir","oakridge1","hubbardbrook","alligatorriver","readingma","bullshoals",
           "willowcreek","downerwoods","laurentides","russellsage","sanford","boundarywaters")
sites <- as.character(siteData$siteName)
avgNum <- 15
cutoff <- -0.02
#n=182
Nmc <- 2000
NT <- 184
days <- seq(1,NT)
offset <- 20
tranOffsets = transData = read.csv('phenocamTransitions_fromMean.csv',header=TRUE)

foreach(calSite = calSites) %dopar% {
  calTran <- as.numeric(tranOffsets[siteData$siteName==calSite,2])
  
  tranID <- which(tranOffsets$siteName==calSite)
  tran_DOY <- as.numeric(tranOffsets[tranID,2])
  
  n=round(tran_DOY+offset,digits=0)-182
  #n=65
  print(n)
  calFileName <- paste0(calSite,"_meanTemp_summer",n,"_expBreak_slope_forecast_b3_calibration_varBurn.RData")
  if(file.exists(calFileName)){
    load(calFileName)
    if(typeof(out.burn)!=typeof(FALSE)){
      out.mat <- data.frame(as.matrix(out.burn$param))
      prow = sample.int(nrow(out.mat),Nmc,replace=TRUE)
      
      b0 <- out.mat$b0[prow]
      b4 <- out.mat$b4[prow]
      b3 <- out.mat$b3[prow]
      b1 <- rep(0,Nmc)
      b2 <- rep(0,Nmc)
      
      outputData <- matrix(nrow=0,ncol=3)
      colnames(outputData) <- c("siteName","year","minDiffDiff_DOY")
      for(s in 1:length(sites)){
        #  for(s in 1:2){
        siteName <- as.character(sites[s])
        print(siteName)
        
        load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData"))
        
        for(yr in 1:dataFinal$N){
          yrName <- dataFinal$years[yr]
          print(yrName)
          initialXs <- rbeta(Nmc,dataFinal$x1.a[yr],dataFinal$x1.b[yr])
          yrPred <- (forecastStep(IC=initialXs,b0=b0,b1=b1,b2=b2,b3=b3,b4=b4,
                                  n=Nmc,NT=NT,Tair=dataFinal$TairMu[,yr],D=dataFinal$D[,yr],Q=sqrt(1/out.mat$p.proc[prow])))
          
          avgDiffs <- matrix(nrow=Nmc,ncol=ncol(yrPred))
          for(g in 1:Nmc){
            avgDiffs[g,] <- calculateStart3(ys=as.numeric(yrPred[g,]),avgNum = avgNum)
          }
          # ci <- apply(avgDiffs,quantile,MARGIN = 2,c(0.025,0.5,0.975))
          # pred.ci <- apply(yrPred,quantile,MARGIN = 2,c(0.025,0.5,0.975))
          # newRow <- c(siteName,yrName,which.min(ci[2,]))
          newRow <- c(siteName,yrName,which.min(colMeans(avgDiffs)))
          outputData <- rbind(outputData,newRow)
        }
      }
      
      save(outputData,file=paste0(calSite,"_inflectionPointData_OOSsites_mean15_",n,".RData"))
    }
  }
}

#dev.off()'
# 
# for(s in 1:16){
#   siteName <- as.character(sites[s])
#   print(siteName)
#   if(file.exists(paste0(siteName,"_inflectionPointData.RData"))){
#   load(paste0(siteName,"_inflectionPointData.RData"))
# 
#   if(s==1){
#     allOutput <- outputData
#   }else{
#     allOutput <- rbind(allOutput,outputData)
#   }
#   }
# }
# 
# #
# # offsetDat <- matrix(nrow=0,ncol=3)
# # colnames(offsetDat) <- c("siteName","year","inflectionOffset")
# outputData$siteYear <- paste0(outputData$siteName,"_",outputData$year)
# 
# outputData <- allOutput
# outputData$siteYear <- paste0(outputData$siteName,"_",outputData$year)
# uniqueSiteYears <- unique(outputData$siteYear)
# reachedCutoffs <- rep(NA,length(uniqueSiteYears))
# reachedCutoffs <- cbind(uniqueSiteYears,reachedCutoffs)
# cutOff <- -0.02
# for(i in 1:length(uniqueSiteYears)){
#   subDat <- outputData[outputData$siteYear==uniqueSiteYears[i],]
#   if(length(which(subDat$minDiffDiff<cutOff))>0){
#   reachedCutoffs[i,2] <- subDat$offset[min(which(subDat$minDiffDiff<cutOff))]
#   }else{
#     reachedCutoffs[i,2] <- NA
#   }
# }
# 
# 
# 
# percentages <- rep(NA,length(offsetVals))
# offsetVals <- seq(min(outputData$offset),max(outputData$offset))
# for(i in 1:length(offsetVals)){
#   percentages[i] <- sum(na.omit(as.numeric(reachedCutoffs[,2]))<=offsetVals[i])/nrow(reachedCutoffs)
# }
# plot(offsetVals,percentages*100,pch=20,typ="l",lwd=2,ylab="Cummulative Percentage of Site-Years (%)",
#      xlab="Included DOYs Relative to Transition DOY (Days)",bty="n")
# 
# OOSdat <- outputData[!as.logical(outputData$includedYear),]
# 
# uniqueSiteYears <- unique(OOSdat$siteYear)
# reachedCutoffs <- rep(NA,length(uniqueSiteYears))
# reachedCutoffs <- cbind(uniqueSiteYears,reachedCutoffs)
# cutOff <- -0.02
# for(i in 1:length(uniqueSiteYears)){
#   subDat <- OOSdat[OOSdat$siteYear==uniqueSiteYears[i],]
#   reachedCutoffs[i,2] <- subDat$offset[min(which(subDat$minDiffDiff<cutOff))]
# }
# 
# percentagesOOS <- rep(NA,length(offsetVals))
# offsetValsOOS <- seq(min(outputData$offset),max(outputData$offset))
# for(i in 1:length(offsetValsOOS)){
#   percentagesOOS[i] <- sum(as.numeric(reachedCutoffs[,2])<=offsetValsOOS[i])/nrow(reachedCutoffs)
# }
# lines(offsetValsOOS,percentagesOOS*100,pch=20,typ="l",lwd=2,col="cyan")
# 
