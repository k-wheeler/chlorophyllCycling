library('rjags')
library('runjags')
library('ecoforecastR')
library('RColorBrewer')
library(doParallel)
source('generalVariables.R')

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
forecastStep <- function(IC,b0,b3,b4,Q=0,n,NT,Tair,D){
  x <- matrix(NA,n,NT)
  if(length(IC)==1){
    Xprev <- rep(IC,n)
  }else{
    Xprev <- IC
  }
  
  for(t in 1:NT){
    bd <- Xprev + b4 * Xprev 
    syn <- b0 + b3 * Tair[t] * D[t] 
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

calculateStart <- function(ys,avgNum=10){
  avgDiffs <- rep(0,avgNum)
  for(t in (avgNum + 1):(length(ys)-avgNum)){
    prevAvg <- mean(diff(ys[1:(t-1)]),na.rm = TRUE)
    newAvg <- mean(diff(ys[t:(t+avgNum)]),na.rm = TRUE)
    avgDiffs <- c(avgDiffs,(newAvg-prevAvg))
  }
  avgDiffs <- c(avgDiffs,rep(0,avgNum))
  return(avgDiffs)
}

registerDoParallel(cores=n.cores)

allSites <- as.character(siteData$siteName)
avgNum <- 15
cutoff <- -0.02

Nmc <- 2000
NT <- 184
days <- seq(1,NT)
offset <- 20
tranOffsets = transData = read.csv(allPhenoTranFile,header=TRUE)

foreach(calSite = sites) %dopar% {
  calTran <- as.numeric(tranOffsets[siteData$siteName==calSite,2])
  
  tranID <- which(tranOffsets$siteName==calSite)
  tran_DOY <- as.numeric(tranOffsets[tranID,2])
  
  n=round(tran_DOY+offset,digits=0)-182
  #n=65
  print(n)
  calFileName <- paste0(CCmodelOutputsFolder,calSite,"_",n,"_ccModel_forecast_calibration_varBurn.RData")
  if(file.exists(calFileName)){
    load(calFileName)
    if(typeof(out.burn)!=typeof(FALSE)){
      out.mat <- data.frame(as.matrix(out.burn$param))
      prow = sample.int(nrow(out.mat),Nmc,replace=TRUE)
      
      b0 <- out.mat$b0[prow]
      b4 <- out.mat$b4[prow]
      b3 <- out.mat$b3[prow]
      
      outputData <- matrix(nrow=0,ncol=3)
      colnames(outputData) <- c("siteName","year","minDiffDiff_DOY")
      for(s in 1:length(allSites)){
        siteName <- as.character(allSites[s])
        print(siteName)
        
        load(paste0(dataDirectory,siteName,"_dataFinal.RData"))
        
        for(yr in 1:dataFinal$N){
          yrName <- dataFinal$years[yr]
          print(yrName)
          initialXs <- rbeta(Nmc,dataFinal$x1.a[yr],dataFinal$x1.b[yr])
          yrPred <- (forecastStep(IC=initialXs,b0=b0,b3=b3,b4=b4,
                                  n=Nmc,NT=NT,Tair=dataFinal$TairMu[,yr],D=dataFinal$D[,yr],Q=sqrt(1/out.mat$p.proc[prow])))
          
          avgDiffs <- matrix(nrow=Nmc,ncol=ncol(yrPred))
          for(g in 1:Nmc){
            avgDiffs[g,] <- calculateStart(ys=as.numeric(yrPred[g,]),avgNum = avgNum)
          }
          newRow <- c(siteName,yrName,which.min(colMeans(avgDiffs)))
          outputData <- rbind(outputData,newRow)
        }
      }
      
      save(outputData,file=paste0(calSite,"_inflectionPointData_OOSsites_mean15_",n,".RData"))
    }
  }
}

