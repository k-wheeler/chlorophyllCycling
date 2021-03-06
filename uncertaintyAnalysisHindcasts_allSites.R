library('rjags')
library('runjags')
library('ecoforecastR')
library('RColorBrewer')
library(scoringRules)

library(doParallel)

n.cores <- 24
registerDoParallel(cores=n.cores)


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


uncertaintyAnalysisHindcast_expBreak <- function(calSite,n=183){
  dataDirectory <- "data/"
  tranOffsets <- read.csv('phenocamTransitions_fromMean.csv',header=TRUE)
  siteData <- read.csv('/projectnb/dietzelab/kiwheel/chlorophyllCycling/allPhenocamDBsitesComplete.csv',header=TRUE)
  Nmc <- 1000
  NT <- 184
  days <- seq(1,NT)
  
  ##Load calibration data:
  calTran <- as.numeric(tranOffsets[siteData$siteName==calSite,2])
  

  tranID <- which(tranOffsets$siteName==calSite)
  tran_DOY <- as.numeric(tranOffsets[tranID,2])
  
  #n=round(tran_DOY+offset,digits=0)-182
  print(n)
  calFileName <- paste0('finalVarBurns/',calSite,"_meanTemp_summer",n,"_expBreak_slope_forecast_b3_calibration_varBurn.RData")
  
  if(file.exists(calFileName)){
    load(calFileName)
    if(typeof(out.burn)==typeof(FALSE)){
      calFileName <- paste0('finalVarBurns/partial_',calSite,"_meanTemp_summer",n,"_expBreak_slope_forecast_b3_calibration_varBurn.RData")
      load(calFileName)
      out.burn <- partialOutput
    }
  }else{
    calFileName <- paste0('finalVarBurns/partial_',calSite,"_meanTemp_summer",n,"_expBreak_slope_forecast_b3_calibration_varBurn.RData")
    load(calFileName)
    out.burn <- partialOutput
  }
  
  out.mat <- data.frame(as.matrix(out.burn$param))
  prow = sample.int(nrow(out.mat),Nmc,replace=TRUE)
  
  b0 <- out.mat$b0[prow]
  b4 <- out.mat$b4[prow]
  b3 <- out.mat$b3[prow]
  b1 <- rep(0,Nmc)
  b2 <- rep(0,Nmc)
  crpsDat <- matrix(ncol=2208+2,nrow=nrow(siteData)*2)
  
  for(s in 1:nrow(siteData)){
    #for(s in 1:3){
    siteName <- as.character(siteData$siteName[s])
    print(siteName)
    if(siteName!="asa2"){
      load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData")) #Load Data
      
      climFileName <- paste0(siteName,"_climatology_forecast_calibration_varBurn2.RData")
      load(climFileName)
      pred.matYr <- data.frame(as.matrix(out.burn))[,1:184]
      pred.mat <- pred.matYr
      if(dataFinal$N>1){
        for(yr in 2:dataFinal$N){
          pred.mat <- cbind(pred.mat,pred.matYr)
        }
      }
      pVals <- as.vector(dataFinal$p)
      crpsClim <- numeric()
      for(i in 1:ncol(pred.mat)){
        if(!is.na(pVals[i])){
          crpsClim <- c(crpsClim,crps_sample(y=pVals[i],dat=pred.mat[,i]))
        }else{
          crpsClim <- c(crpsClim,NA)
        }
      }
      
      initialXs <- rbeta(Nmc,dataFinal$x1.a[1],dataFinal$x1.b[1])
      #print(mean(initialXs))
      ysPred <- forecastStep(IC=initialXs,b0=b0,b1=b1,b2=b2,b3=b3,b4=b4,
                             n=Nmc,NT=NT,Tair=dataFinal$TairMu[,1],D=dataFinal$D[,1],Q=sqrt(1/out.mat$p.proc[prow]))
      if(dataFinal$N>1){
        for(i in 2:dataFinal$N){
          initialXs <- rbeta(Nmc,dataFinal$x1.a[i],dataFinal$x1.b[i])
          ysPred <- cbind(ysPred,forecastStep(IC=initialXs,b0=b0,b1=b1,b2=b2,b3=b3,b4=b4,
                                              n=Nmc,NT=NT,Tair=dataFinal$TairMu[,i],D=dataFinal$D[,i],Q=sqrt(1/out.mat$p.proc[prow])))
        }
      }
      ci = apply(ysPred,2,quantile,c(0.025,0.5,0.975))
      par(mfrow=c(2,1))
      tran_days <- numeric()
      mean_tran_days <- numeric()
      cal_tran_days <- numeric()
      for(yr in 1:dataFinal$N){
        tran_days <- c(tran_days,184*(yr-1) + (as.numeric(tranOffsets[s,(yr+2)]) + as.numeric(tranOffsets[s,2])-182))
        mean_tran_days <- c(mean_tran_days,184*(yr-1) + as.numeric(tranOffsets[s,2])-182)
        cal_tran_days <- c(cal_tran_days,184*(yr-1) + calTran-182)
      }
      
      plot(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,main=paste(siteName,"with calibration",calSite,n),ylim=c(0,1),xlim=c(0,2500))
      ecoforecastR::ciEnvelope(seq(1,length(dataFinal$p)),ci[1,],ci[3,],col=col.alpha("green",0.5))
      abline(v=tran_days,col="red",lty=2,lwd=1)
      abline(v=mean_tran_days,col="red",lwd=1)
      abline(v=cal_tran_days,col="cyan",lwd=1)
      
      #calculate CRPSs
      crpsVals <- numeric()
      pVals <- as.vector(dataFinal$p)
      for(i in 1:ncol(ysPred)){
        if(!is.na(pVals[i])){
          crpsVals <- c(crpsVals,crps_sample(y=pVals[i],dat=ysPred[,i]))
        }else{
          crpsVals <- c(crpsVals,NA)
        }
      }
      crpsDat[(s*2-1),1] <- siteName
      crpsDat[(s*2),1] <- siteName
      crpsDat[(s*2-1),2] <- "crps"
      crpsDat[(s*2),2] <- "crps_diff"
      crpsDat[(s*2-1),3:(length(crpsVals)+2)] <- crpsVals
      crpsDat[(s*2),3:(length(crpsVals)+2)] <- crpsVals-crpsClim
      
      plot(seq(1,length(dataFinal$p)),crpsVals-crpsClim,type="l",main=paste0(siteName," ",n),xlim=c(0,2500))
      abline(h=0,col="gray")
      abline(v=tran_days,col="red",lty=2,lwd=1)
      abline(v=mean_tran_days,col="red",lwd=1)
      abline(v=cal_tran_days,col="cyan",lwd=1)
      #abline(v=seq(n,(dataFinal$n*n),dataFinal$n),col="cyan")
      #lines(seq(1,length(dataFinal$p)),crpsClim,type="l",main=paste0(siteName," ",n),xlim=c(0,2500),col="red")
    }
  }
  write.table(x=crpsDat,sep=",",file=paste0("outOfSampleSites_crps_",calSite,"_183.csv"),row.names=FALSE,col.names=FALSE,quote = FALSE)
}

sites <- c("harvard","umichbiological","bostoncommon","coweeta","howland2",
           "morganmonroe","missouriozarks","queens","dukehw","lacclair","bbc1","NEON.D08.DELA.DP1.00033",
           "bartlettir","oakridge1","hubbardbrook","alligatorriver","readingma","bullshoals",
           "willowcreek","downerwoods","laurentides","russellsage","sanford","boundarywaters")


#output <- read.csv("identifiedDivergenceDOYs.csv",header=TRUE)

foreach(c =1:length(sites)) %dopar% {
#for(c in 1:length(sites)){
  siteName <- sites[c]
  calSite <- siteName
  #if(!file.exists(paste0("outOfSampleSites_crps_",calSite,"_20offset_updated.csv"))){
    pdfName <- paste0("UncertaintyAnaylisHindcast_",calSite,"_Calibration_crps_allSites_183.pdf")
    print(pdfName)
    pdf(file=pdfName,height=10,width=60)
    uncertaintyAnalysisHindcast_expBreak(calSite=calSite)
    dev.off()
  #}
}
# calSite="harvard"
# pdfName <- paste0("UncertaintyAnaylisHindcast_harvardCalibration_91_crps_allSites.pdf")
# pdf(file=pdfName,height=10,width=60)
# n=91
# uncertaintyAnalysisHindcast_expBreak(n=n,calSite=calSite)
# dev.off()
