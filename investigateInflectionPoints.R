library('rjags')
library('runjags')
library('ecoforecastR')
library('RColorBrewer')
library(doParallel)

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

dataDirectory <- "data/"
siteData <- read.csv('/projectnb/dietzelab/kiwheel/chlorophyllCycling/allPhenocamDBsitesComplete.csv',header=TRUE)
sites <- c("harvard","umichbiological","bostoncommon","coweeta","howland2",
           "morganmonroe","missouriozarks","queens","dukehw","lacclair","bbc1","NEON.D08.DELA.DP1.00033",
           "bartlettir","oakridge1","hubbardbrook","alligatorriver","readingma","bullshoals",
           "willowcreek","downerwoods","laurentides","russellsage","sanford","boundarywaters") ##Calibration Sites 

#Not current code, but used to investigate different number of days to average over
# avgNums <- c(10,15,20,25,30,40)
# 
# #pdf(file="investigatingInflectionPoints3.pdf",height=4.5,width=33)
# ##Investigate PhenoCam
# pdf(file="investigatingInflectionPoints.pdf",height=9,width=11)
# par(mfcol=c(2,1))
# par(mai=c(0.5,0.5,0.5,0.1))
# 
# for(s in 1:length(sites)){
# #  for(s in 1:3){
# 
#   siteName <- as.character(sites[s])
#   print(siteName)
# 
#   load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData"))
# 
#   for(yr in 1:dataFinal$N){
#     yrName <- dataFinal$years[yr]
#     # for(n in avgNums){
#     n <- 15
#     avgDiffs <- calculateStart(ys=dataFinal$p[,yr],avgNum = n)
# 
#     #plot(dataFinal$p[,yr],pch=20,main=paste(siteName,yrName,n),ylab="",xlab="")
#     plot(avgDiffs,typ="l",main=paste(siteName,yrName,n),ylab="",xlab="")
#     abline(h=-0.025,col="red")
#     # }
#   }
# }
# dev.off()

##Investigate model predictions!
n.cores <- 16
registerDoParallel(cores=n.cores)
avgNum <- 15
cutoff <- -0.02
#n=182
Nmc <- 1000
#ns <- c(49,56,63,70,77,84,91,98,105,112,119,126,183)
ns <- seq(35,183)

#pdf(file="investigatingInflectionPoints_model.pdf",height=12,width=11)
#par(mfcol=c(3,1))
#par(mai=c(0.5,0.5,0.5,0.1))
transData <- read.csv('phenocamTransitions_fromMean.csv',header=TRUE)
#for(s in 1:length(sites)){
#for(s in 1:16){
foreach(s =1:length(sites)) %dopar% {
  outputData <- matrix(nrow=0,ncol=8)
  colnames(outputData) <- c("siteName","year","includedYear","trans","n","converged","minDiffDiff","meanDiffDiff")
  siteName <- as.character(sites[s])
  print(siteName)
  #if(!file.exists(paste0(siteName,"_inflectionPointData_Updated_10.RData"))){
  
  load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData"))
  siteTrans <- transData[transData$siteName==siteName,]
  for(n in ns){
    fileName <- paste0('finalVarBurns/',siteName,"_meanTemp_summer",n,"_expBreak_slope_forecast_b3_calibration_varBurn.RData")
    partialFileName <- paste0('finalVarBurns/partial_',siteName,"_meanTemp_summer",n,"_expBreak_slope_forecast_b3_calibration_varBurn.RData")
    if(file.exists(fileName)){
      res <- try(load(fileName))#Load Model Output 
      if(inherits(res,"try-error")){
        next
      }
      if(typeof(out.burn)==typeof(FALSE)){
        res <- try(load(partialFileName))#Load Model Output 
        if(inherits(res,"try-error")){
          next
        }
        out.burn <- partialOutput
      }
      
    }else if(file.exists(partialFileName)){
      res <- try(load(partialFileName))#Load Model Output 
      if(inherits(res,"try-error")){
        next
      }
      out.burn <- partialOutput
    }else{
      next 
    }
    pred.mat <- data.frame(as.matrix(out.burn$predict))
    
    prow = sample.int(nrow(pred.mat),Nmc,replace=TRUE)
    pred.mat <- pred.mat[prow,]
    out.mat <- data.frame(as.matrix(out.burn$param))
    
    b0 <- out.mat$b0[prow]
    b4 <- out.mat$b4[prow]
    b3 <- out.mat$b3[prow]
    if(diff(range(b0))<0.5 & diff(range(b4))<0.5 & diff(range(b3))<0.5){
      converged <- TRUE
    }else{
      converged <- FALSE
    }
    #print(converged)
    #print(gelman.diag(partialOutput$params))
    for(yr in 1:dataFinal$N){
      yrName <- dataFinal$years[yr]
      print(yrName)
      if(converged){
        
        yrPred <- pred.mat[1:Nmc,(184*(yr-1)+1):(184*(yr-1)+184)]
        
        avgDiffs <- matrix(nrow=Nmc,ncol=ncol(yrPred))
        for(g in 1:Nmc){
          avgDiffs[g,] <- calculateStart(ys=as.numeric(yrPred[g,]),avgNum = avgNum)
        }
        ci <- apply(avgDiffs,quantile,MARGIN = 2,c(0.025,0.5,0.975))
        pred.ci <- apply(yrPred,quantile,MARGIN = 2,c(0.025,0.5,0.975))
        newRow <- c(siteName,yrName,yrName!=dataFinal$yearRemoved,(siteTrans[yr+2]+siteTrans[2]),
                    n,converged,round(min(ci),digits=4),round(min(ci[2,]),digits=4))
        # plot(dataFinal$p[,yr],pch=20,main=paste(siteName,yrName,n,"PhenoCam"),ylab="",xlab="",ylim=c(0,1))
        # plot(seq(1,ncol(pred.ci)),pred.ci[2,],typ="l",main=paste(siteName,yrName,n,"Predicted"),ylab="",xlab="",ylim=c(0,1))
        # ciEnvelope(seq(1,ncol(pred.ci)),yhi = pred.ci[3,],ylo = pred.ci[1,],col="lightblue")
        # lines(seq(1,ncol(pred.ci)),pred.ci[2,])
        # abline(v=n,col="red")
        # plot(seq(1,ncol(ci)),ci[2,],typ="l",main=paste(siteName,yrName,n),ylab="",xlab="")
        # ciEnvelope(seq(1,ncol(ci)),yhi = ci[3,],ylo = ci[1,],col="lightblue")
        # lines(seq(1,ncol(ci)),ci[2,])
        # abline(v=n,col="red")
        # abline(h=cutoff,col="red")
      }else{
        newRow <- c(siteName,yrName,yrName!=dataFinal$yearRemoved,(siteTrans[yr+2]+siteTrans[2]),n,converged,NA,NA)
      }
      outputData <- rbind(outputData,newRow)
      # }
    }
  }
  outputData <- data.frame(outputData)
  outputData$offset <- round((as.numeric(outputData$n)+182)-as.numeric(outputData$trans),digits=0)
  
  save(outputData,file=paste0(siteName,"_inflectionPointData_15.RData"))
  #}
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
