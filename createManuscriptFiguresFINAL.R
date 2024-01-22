###Paper Figures
#General Library and Data Loading ----
library('rjags')
library('runjags')
library('ecoforecastR')
library('RColorBrewer')
library('scoringRules')
library('scales')
dataDirectory <- "data/"

sites <- c("harvard","umichbiological","bostoncommon","coweeta","howland2",
           "morganmonroe","missouriozarks","queens","dukehw","lacclair","bbc1","NEON.D08.DELA.DP1.00033",
           "bartlettir","oakridge1","hubbardbrook","alligatorriver","readingma","bullshoals",
           "willowcreek","downerwoods","laurentides","russellsage","sanford","boundarywaters") ##Calibration Sites 

siteData <- read.csv('/projectnb/dietzelab/kiwheel/chlorophyllCycling/allPhenocamDBsitesComplete.csv',header=TRUE)

tranOffsets <- read.csv('phenocamTransitions_fromMeanFiltered.csv',header=TRUE)
cutOff <- -0.02 #Cut-off For the value of second-differences that constitutes an inflection based off of PhenoCam data
NT <- 184
avgN <- 15
days <- seq(1,NT)

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
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

##Supplementary Figure to show names of different parts of the curve ----
library(PhenologyBayesModeling)
URL <- "http://phenocam.sr.unh.edu/data/archive/harvard/ROI/harvard_DB_0001_1day.csv"
phenoDataSub <- download.phenocam(URL)

jpeg("curvePartsFigureRAW.jpeg",width=10,height=4.4,units = "in",res=1000)
plot(phenoDataSub$gcc_90[phenoDataSub$year==2013],pch=20,ylab="PhenoCam GCC",xlab="Day of Year",bty="n")
polygon(c(0.01,182,182,0.01),c(0.355,0.355,0.50,0.50),col="gray",border = NA)
points(phenoDataSub$gcc_90[phenoDataSub$year==2013],pch=20)
dev.off()


#Figure 1: Time-series Examples ----
##Example Time-series model fit figure: willowcreek 76; lacclair 43; lacclair 65; howland2 183; alligator 183 
selectSites <- c("howland2","alligatorriver","lacclair","lacclair","willowcreek")
ns <- c(183,183,43,65,76)
tles <- c("Howland, Maine, USA: All Autumn Days/Year Calibration",
          "Alligator River, North Carolina, USA: All Autumn Days/Year Calibration",
          "Lacclair, Quebec, Canada: 43 Days/Year Calibration",
          "Lacclair, Quebec, Canada: 65 Days/Year Calibration",
          "Willow Creek, Wisconsin, USA: 76 Days/Year Calibration")
dataDirectory <- "data/"
Nmc <- 5000

newYears <- numeric()
for(yr in 1:12){
  newYears <- c(newYears,(184*(yr-1)+1))
}
dates <- seq(1,(184*10))
jpeg("timeseries_exampleFigure.jpeg",width=12,height=7.5,units = "in",res=1000)
par(mfrow=c(5,1))
par(mai=c(0.25,0.3,0.5,0.1))
for(i in 1:5){
  siteName <- selectSites[i]
  print(siteName)
  n=ns[i]
  load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData")) #Load Data
  yearRemoved <- dataFinal$yearRemoved
  
  yearInt <- which(dataFinal$years==yearRemoved)
  fileName <- paste0('finalVarBurns/',siteName,"_meanTemp_summer",n,"_expBreak_slope_forecast_b3_calibration_varBurn.RData")
  if(file.exists(fileName)){
    load(fileName)
    if(typeof(out.burn)==typeof(FALSE)){
      fileName <- paste0('finalVarBurns/partial_',siteName,"_meanTemp_summer",n,"_expBreak_slope_forecast_b3_calibration_varBurn.RData")
      load(fileName)
      out.burn <- partialOutput
    }
  }else{
    fileName <- paste0('finalVarBurns/partial_',siteName,"_meanTemp_summer",n,"_expBreak_slope_forecast_b3_calibration_varBurn.RData")
    load(fileName)
    out.burn <- partialOutput
  }

  out.mat <- data.frame(as.matrix(out.burn$params))
  prow = sample.int(nrow(out.mat),Nmc,replace=TRUE)
  
  b0 <- out.mat$b0[prow]
  b4 <- out.mat$b4[prow]
  b3 <- out.mat$b3[prow]
  pred.mat <- data.frame(as.matrix(out.burn$predict))
  ci = apply(pred.mat,2,quantile,c(0.025,0.5,0.975))
  #dates <- as.Date(seq(1,length(dataFinal$p)),origin=as.Date(paste0(dataFinal$years[1],"-01-01")))
  
  if(i==1){
    dataFinal$p <- cbind(dataFinal$p[,1:4],rep(NA,184),dataFinal$p[,5:9])
  }else if(i == 2){
    dataFinal$p <- cbind(rep(NA,184),rep(NA,184),dataFinal$p[,1:5],rep(NA,184),dataFinal$p[,6:7])
  }else if(i == 3 | i==4){
    dataFinal$p <- cbind(rep(NA,184),rep(NA,184),rep(NA,184),rep(NA,184),dataFinal$p)
  }else if(i ==5){
    dataFinal$p <- cbind(rep(NA,184),rep(NA,184),dataFinal$p)
  }
  
  #if(i !=3){
  #  dataFinal$p <- cbind(rep(NA,184),rep(NA,184),dataFinal$p)
  #}
  
  if(i!=5){
    plot(dates,rep(NA,length(dates)),pch=20,xlab="",ylab="",ylim=c(0,1),
         bty="n",main=tles[i],xaxt="n",yaxt="n",xlim=c(0,1850),cex.main=1.5)
    axis(side=1,at=dates[newYears],labels=rep("",12),pos=-0.02,cex=2)
    includeP <- dataFinal$p
    excludeP <- dataFinal$p
    if(n!=183){
      includeP[(n+1):nrow(includeP),] <- NA
    }
    if(i==1){
      excludeP[1:(n+1),-(yearInt+1)] <- NA
      includeP[,(yearInt+1)] <- NA

    }else{ #excluded year comes before any gaps 
      excludeP[1:(n+1),-(yearInt+2)] <- NA
      includeP[,(yearInt+2)] <- NA
    }


    points(dates,includeP,col=col.alpha("black",1),pch=19,cex=0.75)
    
    points(dates,excludeP,col=col.alpha("red",1),pch=19,cex=0.75)
    if(i==1){
      #ecoforecastR::ciEnvelope(dates[369:length(dates)],ci[1,],ci[3,],col=col.alpha("blue",0.50))
      ecoforecastR::ciEnvelope(dates[c(seq(1,4*184),seq((5*184+1),(10*184)))],ci[1,],ci[3,],col=col.alpha("blue",0.50))
      polygon(x=c((4*184),(4*184),(5*184),(5*184)),c(0,1.2,1.2,0),col="white",border=NA)
    }else if(i==2){
      ecoforecastR::ciEnvelope(dates[c(seq((2*184+1),(7*184)),seq((8*184+1),(10*184)))],ci[1,],ci[3,],col=col.alpha("blue",0.50))
      polygon(x=c((7*184),(7*184),(8*184),(8*184)),c(0,1.2,1.2,0),col="white",border=NA)
    }else if(i == 3 | i==4){
      ecoforecastR::ciEnvelope(dates[737:length(dates)],ci[1,],ci[3,],col=col.alpha("blue",0.50))
    }
    #ecoforecastR::ciEnvelope(dates[c(seq(1,4*184),seq((5*184+1),(10*184)))],ci[1,],ci[3,],col=col.alpha("blue",0.50))
  }else{
    plot(dates,rep(NA,length(dates)),pch=20,xlab="Date",ylim=c(0,1),ylab="",
         bty="n",main=tles[i],xaxt="n",yaxt="n",xlim=c(0,1850),cex.main=1.5)
    axis(side=1,at=dates[newYears],labels=seq(as.Date("2011-07-01"),
                                              as.Date("2022-12-31"),"year"),pos=-0.02,cex=2,cex.axis=1.5)
    #points(dates,dataFinal$p,pch=20,col="red")
    includeP <- dataFinal$p
    #includeP[(n+1):nrow(includeP),] <- NA
    includeP[,yearInt] <- NA
    points(dates,includeP,col=col.alpha("black",1),pch=19,cex=0.75)
    
    excludeP <- dataFinal$p
    excludeP[1:n,-yearInt] <- NA
    #excludeP[,-yearInt] <- NA
    points(dates,excludeP,col=col.alpha("red",1),pch=17,cex=0.75)
    ecoforecastR::ciEnvelope(dates[369:length(dates)],ci[1,],ci[3,],col=col.alpha("blue",0.50))
  }
  axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1.0),labels=seq(0,1,0.2),pos=-1,cex=2,cex.axis=1.5)
  if(i==3){
    legend(0,1,col=c("black","red",col.alpha("blue",0.50)),c("Calibration Data","Validation Data","Model Prediction CI"),pch=c(19,17,15),
           bty="n",cex=1.5)
  }
}
dev.off()

#Inflection Figures ----
##Supplementary Figure X: Cumulative ratio of inflection ----

for(s in 1:length(sites)){
  siteName <- as.character(sites[s])
  print(siteName)
  load(paste0(siteName,"_inflectionPointData_15.RData")) #Created in investigateInflectionPoints.R
  
  if(s==1){
    allOutput <- outputData
  }else{
    allOutput <- rbind(allOutput,outputData)
  }
}

outputData <- allOutput
outputData$siteYear <- paste0(outputData$siteName,"_",outputData$year)
uniqueSiteYears <- unique(outputData$siteYear)
reachedCutoffs <- rep(NA,length(uniqueSiteYears))
reachedCutoffs <- cbind(uniqueSiteYears,reachedCutoffs)

for(i in 1:length(uniqueSiteYears)){
  subDat <- outputData[outputData$siteYear==uniqueSiteYears[i],]
  if(length(which(subDat$minDiffDiff<cutOff))>0){
    reachedCutoffs[i,2] <- subDat$offset[min(which(subDat$minDiffDiff<cutOff))]
  }else{
    reachedCutoffs[i,2] <- NA
  }
}

offsetVals <- seq(min(outputData$offset),max(outputData$offset))
percentages <- rep(NA,length(offsetVals))

for(i in 1:length(offsetVals)){
  percentages[i] <- sum(na.omit(as.numeric(reachedCutoffs[,2]))<=offsetVals[i])/nrow(reachedCutoffs)
}
jpeg("cummulativeInflectionFigure.jpeg",width=4.85,height=3.63,units = "in",res=1000)
plot(offsetVals,percentages*100,pch=20,typ="l",lwd=2,ylab="Cummulative Percentage of Site-Years (%)",
     xlab="Days Included Relative to Transition Date",bty="n")
abline(h=percentages[offsetVals==0]*100,lty=3,col="gray")
abline(v=0,lty=3,col="gray")

OOSdat <- outputData[!as.logical(outputData$includedYear),]

uniqueSiteYears <- unique(OOSdat$siteYear)
reachedCutoffs <- rep(NA,length(uniqueSiteYears))
reachedCutoffs <- cbind(uniqueSiteYears,reachedCutoffs)
for(i in 1:length(uniqueSiteYears)){
  subDat <- OOSdat[OOSdat$siteYear==uniqueSiteYears[i],]
  reachedCutoffs[i,2] <- subDat$offset[min(which(subDat$minDiffDiff<cutOff))]
}

percentagesOOS <- rep(NA,length(offsetVals))
offsetValsOOS <- seq(min(outputData$offset),max(outputData$offset))
for(i in 1:length(offsetValsOOS)){
  percentagesOOS[i] <- sum(na.omit(as.numeric(reachedCutoffs[,2]))<=offsetValsOOS[i])/nrow(reachedCutoffs)
}
lines(offsetValsOOS,percentagesOOS*100,pch=20,typ="l",lwd=2,col="red",lty=2)
legend('topleft',col=c("black","red"),c("All","Validation"),lty=c(1,2),
       bty="n",cex=1)
dev.off()


##Figure 2: Correlation between maximum inflections: ----
phenoCamInflectionDat <- (matrix(ncol=3,nrow=0))
colnames(phenoCamInflectionDat) <- c("siteName","year","minDiffDiff_PC") #Calculate minimum second differences for PhenoCam data
for(s in 1:length(sites)){
  siteName <- as.character(sites[s])
  print(siteName)
  load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData"))
  
  for(yr in 1:dataFinal$N){
    yrName <- dataFinal$years[yr]
    avgDiffs <- calculateStart(ys=dataFinal$p[,yr],avgNum = avgN)
    phenoCamInflectionDat <-rbind(phenoCamInflectionDat,c(siteName,yrName,round(min(avgDiffs,na.rm = TRUE),digits=4)))
  }
}
allOutput$siteYear <- paste0(allOutput$siteName,allOutput$year)
phenoCamInflectionDat <- data.frame(phenoCamInflectionDat)
phenoCamInflectionDat$siteYear <- paste0(phenoCamInflectionDat$siteName,phenoCamInflectionDat$year)

offsetVals <- seq(min(outputData$offset),max(outputData$offset))

##Pad allOutput to replace missing offset values at the end with last available data
allOutput <- allOutput[as.logical(allOutput$converged),]
for(s in 1:length(sites)){
  print(sites[s])
  siteDat <- subset(allOutput,siteName==sites[s])
  siteYears <- as.numeric(unique(siteDat$year))
  for(yr in 1:length(siteYears)){
    siteYearDat <- subset(siteDat,year==siteYears[yr])
    for(n in offsetVals){
      if((!(n %in% siteYearDat$offset)) & (n>siteYearDat$offset[nrow(siteYearDat)])){
        newRow <- siteYearDat[nrow(siteYearDat),]
        newRow[5] <- as.numeric(newRow[5][[1]]) + (n- as.numeric(newRow[9][[1]]))
        newRow[9] <- n
        allOutput <- rbind(allOutput,newRow)
      }else if(!(n %in% siteYearDat$offset)){
        foundPrevious <- FALSE
        newN <- n
        while(!foundPrevious & newN>=siteYearDat$offset[1]){
          newN <- newN - 1
          if(newN %in% siteYearDat$offset){
            foundPrevious <- TRUE
            newRow <- siteYearDat[siteYearDat$offset==newN,]
            newRow[5] <- as.numeric(newRow[5][[1]]) + (n- as.numeric(newRow[9][[1]]))
            newRow[9] <- n
            allOutput <- rbind(allOutput,newRow)
          }
        }
      }
    }
  }
}

ps <- numeric()
psOOS <- numeric()
fs <- numeric()
fsOOS <- numeric()
sizes <- numeric()
sizesOOS <- numeric()
for(n in offsetVals){
  print(n)
  allOutputTrans <- allOutput[allOutput$offset==n,] 
  mergedDat <- merge(allOutputTrans,phenoCamInflectionDat)
  mergedDat$minDiffDiff_PC <- as.numeric(as.character(mergedDat$minDiffDiff_PC))
  if(length(na.omit(as.numeric(mergedDat$minDiffDiff)))>1){
    mdl <- lm(as.numeric(mergedDat$minDiffDiff_PC)~as.numeric(mergedDat$minDiffDiff))
    sm <- summary(mdl)
    ps <- c(ps,lmp(mdl))
    fs <- c(fs,sm$fstatistic[1])
    sizes <- c(sizes,length(na.omit(as.numeric(mergedDat$minDiffDiff))))
  }else{
    ps <- c(ps,NA)
    fs <- c(fs,NA)
    sizes <- c(sizes,NA)
  }
  validationDat <- mergedDat[!as.logical(mergedDat$includedYear),]
  if(length(na.omit(as.numeric(validationDat$minDiffDiff)))>2){
    mdl <- lm(as.numeric(validationDat$minDiffDiff_PC)~as.numeric(validationDat$minDiffDiff))
    sm <- summary(mdl)
    psOOS <- c(psOOS,lmp(mdl))
    fsOOS <- c(fsOOS,sm$fstatistic[1])
    sizesOOS <- c(sizesOOS,length(na.omit(as.numeric(validationDat$minDiffDiff))))
  }else{
    psOOS <- c(psOOS,NA)
    fsOOS <- c(fsOOS,NA)
    sizesOOS <- c(sizesOOS,NA)
    
  }
}
jpeg("correlationPvsIncludedDOY_Figure_andExample.jpeg",width=8,height=4,units = "in",res=1000)
par(mai=c(0.8,0.8,0.5,0.8))
par(mfrow=c(1,2))

n=10
#n=127
print(n)
allOutputTrans <- allOutput[allOutput$offset==n,]
mergedDat <- merge(allOutputTrans,phenoCamInflectionDat)
mergedDat$minDiffDiff_PC <- as.numeric(as.character(mergedDat$minDiffDiff_PC))

mdl <- lm(as.numeric(mergedDat$minDiffDiff_PC)~as.numeric(mergedDat$minDiffDiff))
sm <- summary(mdl)
plot(as.numeric(mergedDat$minDiffDiff),as.numeric((mergedDat$minDiffDiff_PC)),pch=20,ylim=c(-0.08,0.01),xlab="Model Second Difference",
     ylab="PhenoCam Second Difference",bty="n")
abline(sm$coefficients[1],sm$coefficients[2],col="cyan",lwd=3,lty=4)
legend("topleft",lwd=3,lty=4,"Line of Best Fit",col="cyan",cex=1.25,bty="n")

plot(offsetVals,ps,pch=20,bty="n",xlab="Days Included Relative to Transition Date",ylab="p-value",xlim=c(-70,130),ylim=c(0,1))
abline(h=0.05,col="gray",lty=2,lwd=2)
points(offsetVals,ps,pch=20)
#points(offsetVals,psOOS,pch=20,col=alpha("red",0.5))
points(10,ps[offsetVals==10],col="cyan",pch=20)
#points(offsetVals,psOOS,col="red",pch=17)
# legend('topleft',col=c("black","red"),c("All","Validation"),pch=c(19,17),
#        bty="n",cex=1)
dev.off()

##Supplementary Figure **: Statistics and include validation years 
jpeg("correlationPvsIncludedDOY_statisticsWithValidation.jpeg",width=5,height=5,units = "in",res=1000)
par(mai=c(0.6,0.6,0.2,0.2))
par(mfrow=c(3,1))
plot(offsetVals,ps,pch=20,bty="n",xlab="",ylab="p-value",xlim=c(-50,100),ylim=c(0,1))
abline(h=0.05,col="gray",lty=2,lwd=2)
points(offsetVals,ps,pch=20)
points(offsetVals,psOOS,pch=17,col=alpha("red",1))
legend('topright',c("All Site-Years","Withheld Years"),pch=c(20,17),col=c("black","red"),bty="n")
text(-50,0.8,"a)",cex=2)

plot(offsetVals,fs,pch=20,bty="n",xlab="",ylab="F-statistic",xlim=c(-50,100))
points(offsetVals,fsOOS,pch=17,col=alpha("red",1))
text(-50,20,"b)",cex=2)

plot(offsetVals,sizes,pch=20,bty="n",xlab="Calibration Days Relative to SOS Date",
     ylab="Sample Size",xlim=c(-50,100))
points(offsetVals,sizesOOS,pch=17,col=alpha("red",1))
text(-50,150,"c)",cex=2)

dev.off()


#Figure 3: CRPS Over time ----

###Need to create 5 rounds of this: 183 full; -8 offset full; 183 sections; -8 sections; validation sections  
hLine=0
ylimVals <- c(-0.2,0.4)
for(rnd in 3:5){
  print(rnd)
  if(rnd==1){ #Start file for main-text figure
    jpeg("hindcastsDifferentWindows_forecasts_crps_diffMain.jpeg",width=6.5,height=9,units = "in",res=1000)
    par(mai=c(0.5,0.8,0.5,0.2))
    par(mfrow=c(2,1))
  }else if(rnd==3){
    jpeg("hindcastsDifferentWindows_forecasts_crps_full.jpeg",width=20,height=12,units = "in",res=1000)
    par(mai=c(0.5,0.8,0.5,0.2))
    par(mfrow=c(3,3))
  }
  if(rnd %in% c(1,3,5)){
    n=183#########*************
    offset <- 183
  }else{
    offset <- -8
    n=round(tran_DOY+offset,digits=0)-182
  }
  allCRPSdiffs <- matrix(nrow=0,ncol=5)
  colnames(allCRPSdiffs) <- c("siteName","year","n_offset","day_offset","val")
  
  crpsDiffsALL <- matrix(ncol=186,nrow=0)
  daysALL <- matrix(ncol=186,nrow=0)
  allTrans <- numeric()
  for(s in 1:length(sites)){
  #for(s in 1:4){
    siteName <- sites[s]
    
    print(siteName)
    load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData")) #Load Data
    load(paste0(siteName,"_inflectionPointData_15.RData"))#outputData
    minDiffDiffDat <- subset(outputData,offset==offset)#outputData[outputData$offset==offset,]
    yearRemoved <- dataFinal$yearRemoved
    
    yearInt <- which(dataFinal$years==yearRemoved)
    climFileName <- paste0(siteName,"_climatology_forecast_calibration_varBurn2.RData")
    load(climFileName)
    pred.matYr <- data.frame(as.matrix(out.burn))[,1:184]
    pred.mat <- pred.matYr
    for(yr in 2:dataFinal$N){
      pred.mat <- cbind(pred.mat,pred.matYr)
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
    if(rnd==5){
      years=dataFinal$yearRemoved
    }else{
      years=dataFinal$years
    }
    
    for(yr in 1:length(years)){
      #print(dataFinal$years[yr])
      yrName=dataFinal$years[yr]
      print(yrName)
      tranID <- which(tranOffsets$siteName==siteName)
      tran_DOY <- as.numeric(tranOffsets[tranID,(yr+2)]) + as.numeric(tranOffsets[tranID,2])
      if(sum(minDiffDiffDat$year==yrName)>0){
        minDiffDiff <- minDiffDiffDat$minDiffDiff[minDiffDiffDat$year==yrName]
        if(is.na(minDiffDiff)){
          minDiffDiff <- 0
        }
      }else{
        minDiffDiff <- 0
      }
      fileName <- paste0('finalVarBurns/',siteName,"_meanTemp_summer",n,"_expBreak_slope_forecast_b3_calibration_varBurn.RData")
      #print(fileName)
      if(file.exists(fileName)){
        load(fileName) #Load Model Output 
        if(typeof(out.burn)!=typeof(FALSE)){
          out.mat <- data.frame(as.matrix(out.burn$param))
          b0 <- out.mat$b0
          b4 <- out.mat$b4
          b3 <- out.mat$b3
          
          if(diff(quantile(b0,c(0.025,0.975)))<0.25 & diff(quantile(b4,c(0.025,0.975)))<0.25 & diff(quantile(b3,c(0.025,0.975)))<0.25){
            convergedWell <- TRUE
          }else{
            convergedWell <- FALSE
          }
          
          pred.mat <- data.frame(as.matrix(out.burn$predict))
          pred.mat.year <- pred.mat[,(184*(yr-1)+1):(184*(yr-1)+184)]
          
          crpsVals <- numeric()
          for(i in 1:ncol(pred.mat.year)){
            if(!is.na(dataFinal$p[i,yr])){
              crpsVals <- c(crpsVals,crps_sample(y=dataFinal$p[i,yr],dat=pred.mat.year[,i]))
            }else{
              crpsVals <- c(crpsVals,NA)
            }
          }
          crps_diff <- crpsVals-(crpsClim[(184*(yr-1)+1):(184*(yr-1)+184)])
          if(convergedWell){
            crpsDiffsALL <- rbind(crpsDiffsALL,c(siteName,yrName,crps_diff))
            daysALL <- rbind(daysALL,c(siteName,yrName,((182 +days)-tran_DOY)))
            allTrans <- c(allTrans,as.numeric(tranOffsets[tranID,(yr+2)]))
          }
        }else{
          fileName <- paste0('finalVarBurns/partial_',siteName,"_meanTemp_summer",n,"_expBreak_slope_forecast_b3_calibration_varBurn.RData")
          if(file.exists(fileName)){
            load(fileName) #Load Model Output 
            out.mat <- data.frame(as.matrix(partialOutput$param))
            b0 <- out.mat$b0
            b4 <- out.mat$b4
            b3 <- out.mat$b3
            
            if(diff(quantile(b0,c(0.025,0.975)))<0.25 & diff(quantile(b4,c(0.025,0.975)))<0.25 & diff(quantile(b3,c(0.025,0.975)))<0.25){
              convergedWell <- TRUE
            }else{
              convergedWell <- FALSE
            }
            
            pred.mat <- data.frame(as.matrix(partialOutput$predict))
            pred.mat.year <- pred.mat[,(184*(yr-1)+1):(184*(yr-1)+184)]
            
            crpsVals <- numeric()
            for(i in 1:ncol(pred.mat.year)){
              if(!is.na(dataFinal$p[i,yr])){
                crpsVals <- c(crpsVals,crps_sample(y=dataFinal$p[i,yr],dat=pred.mat.year[,i]))
              }else{
                crpsVals <- c(crpsVals,NA)
              }
            }
            crps_diff <- crpsVals-(crpsClim[(184*(yr-1)+1):(184*(yr-1)+184)])
            if(convergedWell){
              crpsDiffsALL <- rbind(crpsDiffsALL,c(siteName,yrName,crps_diff))
              daysALL <- rbind(daysALL,c(siteName,yrName,((182 +days)-tran_DOY)))
              allTrans <- c(allTrans,as.numeric(tranOffsets[tranID,(yr+2)]))
            }
          }
        }
      }else{
        fileName <- paste0('finalVarBurns/partial_',siteName,"_meanTemp_summer",n,"_expBreak_slope_forecast_b3_calibration_varBurn.RData")
        if(file.exists(fileName)){
          res <- try(load(fileName))#Load Model Output 
          if(inherits(res,"try-error")){
            next
          }#Load Model Output 
          out.mat <- data.frame(as.matrix(partialOutput$param))
          b0 <- out.mat$b0
          b4 <- out.mat$b4
          b3 <- out.mat$b3
          
          if(diff(quantile(b0,c(0.025,0.975)))<0.25 & diff(quantile(b4,c(0.025,0.975)))<0.25 & diff(quantile(b3,c(0.025,0.975)))<0.25){
            convergedWell <- TRUE
          }else{
            convergedWell <- FALSE
          }
          
          pred.mat <- data.frame(as.matrix(partialOutput$predict))
          pred.mat.year <- pred.mat[,(184*(yr-1)+1):(184*(yr-1)+184)]
          
          crpsVals <- numeric()
          for(i in 1:ncol(pred.mat.year)){
            if(!is.na(dataFinal$p[i,yr])){
              crpsVals <- c(crpsVals,crps_sample(y=dataFinal$p[i,yr],dat=pred.mat.year[,i]))
            }else{
              crpsVals <- c(crpsVals,NA)
            }
          }
          crps_diff <- crpsVals-(crpsClim[(184*(yr-1)+1):(184*(yr-1)+184)])
          if(convergedWell){
            crpsDiffsALL <- rbind(crpsDiffsALL,c(siteName,yrName,crps_diff))
            daysALL <- rbind(daysALL,c(siteName,yrName,((182 +days)-tran_DOY)))
            allTrans <- c(allTrans,as.numeric(tranOffsets[tranID,(yr+2)]))
          }
        }
      }
    }
  }
  
  daysOffset <- seq(round(min((as.numeric(daysALL[,3:ncol(daysALL)])),na.rm=TRUE),digits=0),
                    round(max((as.numeric(daysALL[,3:ncol(daysALL)])),na.rm=TRUE),digits=0))
  reorganizedDat <- data.frame(matrix(nrow=nrow(daysALL),ncol=2+length(daysOffset)))
  reorganizedDat[,1] <- daysALL[,1]
  reorganizedDat[,2] <- daysALL[,2]
  
  for(i in 1:nrow(reorganizedDat)){
    #print(i)
    for(d in 1:length(daysOffset)){
      ind <- which(round(as.numeric(daysALL[i,]),digits=0)==daysOffset[d])
      if(length(ind)>0){
        ind <- ind[1]
        reorganizedDat[i,(d+2)] <- as.numeric(crpsDiffsALL[i,ind])
      }
    }
  }
  
  ps <- numeric()
  mns <- numeric()
  ciTop <- numeric()
  ciBot <- numeric()
  for(j in 3:(ncol(reorganizedDat))){
    if(length(na.omit(reorganizedDat[,j]))>2){
      mdl <- (t.test(na.omit(reorganizedDat[,j])))
      ps <- c(ps,mdl$p.value)
      mns <- c(mns,mdl$estimate)
      ciTop <- c(ciTop,mdl$conf.int[2])
      ciBot <- c(ciBot,mdl$conf.int[1])
    }else{
      ps <- c(ps,NA)
      ciTop <- c(ciTop,NA)
      ciBot <- c(ciBot,NA)
      mns <- c(mns,NA)
    }
  }
  
  plot(numeric(),numeric(),ylim=ylimVals,xlim=c(-125,100),cex.axis=2,bty="n",ylab="",
       main="All Site-Years",cex.main=2,xlab="")
  ecoforecastR::ciEnvelope(daysOffset[!is.na(ciTop)],as.numeric(na.omit(ciBot)),as.numeric(na.omit(ciTop)),col="gray")
  abline(h=hLine,col="black",lty=2,lwd=2)
  lines(daysOffset,mns,col="blue",lwd=3)
  
  if(rnd %in% 3:5){
    earlyDat <- reorganizedDat[allTrans<(-3),]
    ps <- numeric()
    mns <- numeric()
    ciTop <- numeric()
    ciBot <- numeric()
    for(j in 3:(ncol(earlyDat))){
      if(length(na.omit(earlyDat[,j]))>1){
        mdl <- (t.test(na.omit(earlyDat[,j])))
        ps <- c(ps,mdl$p.value)
        mns <- c(mns,mdl$estimate)
        ciTop <- c(ciTop,mdl$conf.int[2])
        ciBot <- c(ciBot,mdl$conf.int[1])
      }else{
        ps <- c(ps,NA)
        ciTop <- c(ciTop,NA)
        ciBot <- c(ciBot,NA)
        mns <- c(mns,NA)
      }
    }
    
    plot(numeric(),numeric(),ylim=ylimVals,xlim=c(-125,100),cex.axis=2,bty="n",ylab="",
         main="Transition Earlier Than Average",cex.main=2,xlab="")
    ecoforecastR::ciEnvelope(daysOffset[!is.na(ciTop)],as.numeric(na.omit(ciBot)),as.numeric(na.omit(ciTop)),col="gray")
    abline(h=hLine,col="black",lty=2,lwd=2)
    lines(daysOffset,mns,col="blue",lwd=3)
    
    
    lateDat <- reorganizedDat[allTrans>3,]
    ps <- numeric()
    mns <- numeric()
    ciTop <- numeric()
    ciBot <- numeric()
    for(j in 3:(ncol(lateDat))){
      if(length(na.omit(lateDat[,j]))>1){
        mdl <- (t.test(na.omit(lateDat[,j])))
        ps <- c(ps,mdl$p.value)
        mns <- c(mns,mdl$estimate)
        ciTop <- c(ciTop,mdl$conf.int[2])
        ciBot <- c(ciBot,mdl$conf.int[1])
      }else{
        ps <- c(ps,NA)
        ciTop <- c(ciTop,NA)
        ciBot <- c(ciBot,NA)
        mns <- c(mns,NA)
      }
    }
    
    plot(numeric(),numeric(),ylim=ylimVals,xlim=c(-125,100),cex.axis=2,bty="n",ylab="",
         main="Transition Later Than Average",cex.main=2,xlab="")
    ecoforecastR::ciEnvelope(daysOffset[!is.na(ciTop)],as.numeric(na.omit(ciBot)),as.numeric(na.omit(ciTop)),col="gray")
    abline(h=hLine,col="black",lty=2,lwd=2)
    lines(daysOffset,mns,col="blue",lwd=3)
  }
  if(rnd %in% c(2,5)){
    dev.off()
  }
}
##Heat map (Figure 3c)
jpeg("CRPS_heatmap_includedVsDay.jpeg",width=4.5,height=4.5,units = "in",res=1000)
par(mfrow=c(1,1))
par(mai=c(0.8,0.8,0.2,0.4))
par(pty="s")
load("crpsMat_includedVsDay_Complete.RData") #loaded as crpsMat; Created in createCRPSpercentageMatrix.R
load(file="daysOffset_includedVsDay.RData") #loaded as daysOffset; Created in createCRPSpercentageMatrix.R
load(file="reorganizedDat_includedVsDay.RData")#loaded as reorganizedDat; Created in createCRPSpercentageMatrix.R
image(x=daysOffset,y=daysOffset,z=t(crpsMat),col=c("#a6611a","#dfc27d","#80cdc1","#018571"),xlab="Days Included Relative to Transition Date",
      ylab="Predicted Day Relative to Transition Date",main="")
segments(0,daysOffset[1],0,daysOffset[length(daysOffset)],col="black",lwd=2,lty=2)
segments(daysOffset[1],0,daysOffset[length(daysOffset)],0,col="black",lwd=2,lty=2)  
segments(daysOffset[1],daysOffset[1],
         daysOffset[length(daysOffset)],daysOffset[length(daysOffset)],col="black",lwd=2,lty=2)
#segments(-8,daysOffset[1],-8,daysOffset[length(daysOffset)],col="red",lwd=3,lty=1)
segments(max(daysOffset),daysOffset[1],max(daysOffset),daysOffset[length(daysOffset)],col="red",lwd=3,lty=1)
legend('topleft',c('< 25%','25-49%','50-74%','> 74%'),
       col=c("#a6611a","#dfc27d","#80cdc1","#018571"),pch=rep(15,4),bty="n")

dev.off()

#Figure 4: Transferability ---- 
plotOOSWithheldYear <- function(type,vl,sites,allSites){
  OOScrps <- matrix(nrow=length(sites),ncol=70)
  temps <- matrix(nrow=length(sites),ncol=70)
  precips <- matrix(nrow=length(sites),ncol=70)
  lats <- matrix(nrow=length(sites),ncol=70)
  mins <- matrix(nrow=length(sites),ncol=70)
  maxs <- matrix(nrow=length(sites),ncol=70)
  dtes <- matrix(nrow=length(sites),ncol=70)
  elevs <- matrix(nrow=length(sites),ncol=70)
  pValsMat <- matrix(nrow=length(sites),ncol=7)
  colnames(pValsMat) <- c("temps","precips","lats","mins","maxs","dtes","elevs")
  tranOffsets <- read.csv('phenocamTransitions_fromMean.csv',header=TRUE)
  
  for(i in 1:length(sites)){

    siteName <- sites[i]
    print(siteName)
    calSite <- siteName
    crpsDat <- read.csv(file=paste0("outOfSampleSites_crps_",calSite,"_183.csv"),header=FALSE) #Created in uncertaintyAnalysisHindcasts_allSites.R

    j <- 1
    calInfo <- siteData[as.character(siteData$siteName)==calSite,]
    
    for(st in 1:length(allSites)){
      siteName <- allSites[st]
      #print(siteName)
      s <- which(crpsDat[,1]==siteName)[vl]
      siteInfo <- siteData[as.character(siteData$siteName)==siteName,]
      load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData")) #Load Data
      yearRemoved <- dataFinal$yearRemoved
      
      yr <- which(dataFinal$years==yearRemoved)
      tranID <- which(tranOffsets$siteName==siteName)
      tran_DOY <- round((as.numeric(tranOffsets[tranID,(yr+2)])+as.numeric(tranOffsets[tranID,2])-182),digits=0)
      
      siteYearDat <- crpsDat[s,((184*(yr-1)+1)+2):((184*(yr-1)+184)+2)]

      if(type=="all"){
        if(vl==1){ #crps vals
          OOScrps[i,j] <- as.numeric((mean(as.numeric(siteYearDat),na.rm=TRUE)))
        }else if(vl==2){ #crps differences
          OOScrps[i,j] <- as.numeric((mean(as.numeric(siteYearDat),na.rm=TRUE)<0))
        }
      }else if(type=="divergence"){
        if(vl==1){ #crps vals
          OOScrps[i,j] <- as.numeric((mean(as.numeric(siteYearDat[1:tran_DOY+20]),na.rm=TRUE)))
        }else if(vl==2){ #crps differences
          OOScrps[i,j] <- as.numeric((mean(as.numeric(siteYearDat[1:tran_DOY+20]),na.rm=TRUE)<0))
        }
      }else if(type=="transition"){
        if(vl==1){ #crps vals
          OOScrps[i,j] <- as.numeric((mean(as.numeric(siteYearDat[(tran_DOY-3):(tran_DOY+3)]),na.rm=TRUE)))
        }else if(vl==2){ #crps differences
          OOScrps[i,j] <- as.numeric((mean(as.numeric(siteYearDat[(tran_DOY-3):(tran_DOY+3)]),na.rm=TRUE)<0))
        }
      }else if(type=="before"){
        if(vl==1){ #crps vals
          OOScrps[i,j] <- as.numeric((mean(as.numeric(siteYearDat[(tran_DOY-6):(tran_DOY)]),na.rm=TRUE)))
        }else if(vl==2){ #crps differences
          OOScrps[i,j] <- as.numeric((mean(as.numeric(siteYearDat[(tran_DOY-6):(tran_DOY)]),na.rm=TRUE)<0))
        }
      }else if(type=="after"){
        if(vl==1){ #crps vals
          OOScrps[i,j] <- as.numeric((mean(as.numeric(siteYearDat[(tran_DOY):(tran_DOY+6)]),na.rm=TRUE)))
        }else if(vl==2){ #crps differences
          OOScrps[i,j] <- as.numeric((mean(as.numeric(siteYearDat[(tran_DOY):(tran_DOY+6)]),na.rm=TRUE)<0))
        }
      }
      temps[i,j] <- abs(siteInfo$MAT - calInfo$MAT)
      precips[i,j] <- abs(siteInfo$MAP - calInfo$MAP)
      lats[i,j] <- abs(siteInfo$Lat - calInfo$Lat)
      mins[i,j] <- abs(siteInfo$rescalingMin - calInfo$rescalingMin)
      maxs[i,j] <- abs(siteInfo$rescalingMax - calInfo$rescalingMax)
      dtes[i,j] <- abs(siteInfo$avgTran - calInfo$avgTran)
      elevs[i,j] <- abs(siteInfo$elevation - calInfo$elevation)
      
      j=j+1
    }
    # if(type=="after"){
    #   print("entered")
    #   #mdl <- lm(OOScrps[i,]~temps[i,]+precips[i,]+lats[i,]+mins[i,]+maxs[i,]+dtes[i,]+elevs[i,])
    #   mdl <- lm(OOScrps[i,]~temps[i,])
    #   sm <- summary(mdl) 
    #   pValsMat[i,1] <- round(sm$coefficients[2,4],digits=3)
    #   mdl <- lm(OOScrps[i,]~precips[i,])
    #   sm <- summary(mdl) 
    #   pValsMat[i,2] <- round(sm$coefficients[2,4],digits=3)
    #   mdl <- lm(OOScrps[i,]~lats[i,])
    #   sm <- summary(mdl) 
    #   pValsMat[i,3] <- round(sm$coefficients[2,4],digits=3)
    #   mdl <- lm(OOScrps[i,]~mins[i,])
    #   sm <- summary(mdl) 
    #   pValsMat[i,4] <- round(sm$coefficients[2,4],digits=3)
    #   mdl <- lm(OOScrps[i,]~maxs[i,])
    #   sm <- summary(mdl) 
    #   pValsMat[i,5] <- round(sm$coefficients[2,4],digits=3)
    #   mdl <- lm(OOScrps[i,]~dtes[i,])
    #   sm <- summary(mdl) 
    #   pValsMat[i,6] <- round(sm$coefficients[2,4],digits=3)
    #   mdl <- lm(OOScrps[i,]~elevs[i,])
    #   sm <- summary(mdl) 
    #   pValsMat[i,7] <- round(sm$coefficients[2,4],digits=3)
    # }
  }
  if(type=="after"){
      #mdl <- lm(OOScrps~temps+precips+lats+mins+maxs[i,]+dtes[i,]+elevs[i,])
        #mdl <- lm(OOScrps[i,]~temps[i,]+precips[i,]+lats[i,]+mins[i,]+maxs[i,]+dtes[i,]+elevs[i,])
        # mdl <- lm(as.vector(OOScrps)~as.vector(temps))
        # sm <- summary(mdl)
        # print("temps")
        # sm$adj.r.squared
        # sm$coefficients[2,4]
        # 
        # mdl <- lm(as.vector(OOScrps)~as.vector(precips))
        # sm <- summary(mdl)
        # print("precips")
        # sm$adj.r.squared
        # sm$coefficients[2,4]
        # 
        # mdl <- lm(as.vector(OOScrps)~as.vector(lats))
        # sm <- summary(mdl)
        # print("lats")
        # sm$adj.r.squared
        # sm$coefficients[2,4]
        # 
        # mdl <- lm(as.vector(OOScrps)~as.vector(mins))
        # sm <- summary(mdl)
        # print("mins")
        # sm$adj.r.squared
        # sm$coefficients[2,4]
        # 
        # mdl <- lm(as.vector(OOScrps)~as.vector(maxs))
        # sm <- summary(mdl)
        # print("maxs")
        # sm$adj.r.squared
        # sm$coefficients[2,4]
        # 
        # mdl <- lm(as.vector(OOScrps)~as.vector(dtes))
        # sm <- summary(mdl)
        # print("dtes")
        # sm$adj.r.squared
        # sm$coefficients[2,4]
        # 
        # mdl <- lm(as.vector(OOScrps)~as.vector(elevs))
        # sm <- summary(mdl)
        # print("elevs")
        # sm$adj.r.squared
        # sm$coefficients[2,4]
    
    # mdl <- lm(as.vector(OOScrps)~as.vector(temps)+as.vector(lats))
    # sm <- summary(mdl)
    # print("temps + lats")
    # sm$adj.r.squared
    # sm$coefficients[,4]
    # 
    # mdl <- lm(as.vector(OOScrps)~as.vector(temps)+as.vector(maxs))
    # sm <- summary(mdl)
    # print("temps + maxs")
    # sm$adj.r.squared
    # sm$coefficients[,4]
    
    mdl <- lm(as.vector(OOScrps)~as.vector(temps)+as.vector(maxs)+as.vector(lats))
    sm <- summary(mdl)
    print("temps + maxs + lats")
    sm$adj.r.squared
    sm$coefficients[,4]
    
  }
  return(OOScrps)
}
OOScrpsALL <- plotOOSWithheldYear(type="all",vl=2,sites=sites,allSites = allSites)
OOScrpsTRAN <- plotOOSWithheldYear(type="transition",vl=2,sites=sites,allSites = allSites)
OOScrpsBEF <- plotOOSWithheldYear(type="before",vl=2,sites=sites,allSites = allSites)
OOScrpsAFT <- plotOOSWithheldYear(type="after",vl=2,sites=sites,allSites = allSites)

cols <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e")
jpeg("transferabilityCount_examplePiePlots.jpeg",width=4.5,height=4.5,units = "in",res=1000)
par(mfrow=c(2,2))
par(mai=c(0.1,0.3,0.3,0.3))
selectedSite <- "umichbiological"
slices <- c(sum(OOScrpsALL[sites==selectedSite ,]),70-sum(OOScrpsALL[sites==selectedSite ,]))
pie(slices,labels=c("% Better","% Worst"),col=c(cols[1],"white"),main="a) Full Autumn")

slices <- c(sum(OOScrpsTRAN[sites==selectedSite ,]),70-sum(OOScrpsTRAN[sites==selectedSite ,]))
pie(slices,labels=c("% Better","% Worst"),col=c(cols[3],"white"),main="b) Within 3 Days")

slices <- c(sum(OOScrpsBEF[sites==selectedSite ,]),70-sum(OOScrpsBEF[sites==selectedSite ,]))
pie(slices,labels=c("% Better","% Worst"),col=c(cols[4],"white"),main="c) 0-6 Days Before")

slices <- c(sum(OOScrpsAFT[sites==selectedSite ,]),70-sum(OOScrpsAFT[sites==selectedSite ,]))
pie(slices,labels=c("% Better","% Worst"),col=c(cols[5],"white"),main="d) 0-6 Days After")

dev.off()

jpeg("transferabilityCount_DensityPlots.jpeg",width=6,height=5,units = "in",res=1000)
plot(density(apply(OOScrpsALL,MARGIN = 1,FUN=sum)/70*100),bty="n",xlab="Percentage of Sites (%)",ylab="Density",
     main="",lwd=2,xlim=c(0,100),col=col.alpha(cols[1],0.4),lty=3)
lines(density(apply(OOScrpsTRAN,MARGIN = 1,FUN=sum,na.rm=TRUE)/70*100),lty=4,lwd=2,col=cols[3])
lines(density(apply(OOScrpsBEF,MARGIN = 1,FUN=sum,na.rm=TRUE)/70*100),lty=5,lwd=2,col=cols[4])
lines(density(apply(OOScrpsAFT,MARGIN = 1,FUN=sum)/70*100),lty=1,lwd=2,col=cols[5])
legend("topright",c("Full Autumn","Within 3 Days","6 Days Before","6 Days After"),
       lwd=rep(2,2),lty=c(3,4,5,1),bty="n",col=cols[c(1,3,4,5)])
dev.off()

output <- cbind(sites,round(apply(OOScrpsAFT,MARGIN = 1,FUN=sum)/70*100,digits=0),
                round(apply(OOScrpsTRAN,MARGIN = 1,FUN=sum,na.rm=TRUE)/70*100,digits=0),
                round(apply(OOScrpsBEF,MARGIN = 1,FUN=sum,na.rm=TRUE)/70*100,digits=0),
                round(apply(OOScrpsALL,MARGIN = 1,FUN=sum)/70*100,digits=0))
output <- output[order(as.numeric(output[,2]),decreasing = TRUE),]
write.csv(output,file="transferabilityPercentBetter.csv",row.names=FALSE,quote=FALSE)

#Supplementary Figure: Boston Common Time-series Examples and Fits ----
library(PhenologyBayesModeling)
siteName <- "bostoncommon"
load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData")) #Load Data

#jpeg("curvePartsFigureRAW.jpeg",width=10,height=4.4,units = "in",res=1000)
plot(as.vector(dataFinal$p),pch=20)
plot(phenoDataSub$date,phenoDataSub$gcc_90,pch=20,ylab="PhenoCam GCC",bty="n")
polygon(c(0.01,182,182,0.01),c(0.355,0.355,0.50,0.50),col="gray",border = NA)
points(phenoDataSub$gcc_90[phenoDataSub$year==2013],pch=20)
dev.off()
