#Investigate Importance of Timing of Environmental Covariates on Predicting Transition
library(scales)
library('mgcv')
library(tidyverse)
library(randomForest)
library(caTools)
library(tidyverse)
source('generalVariables.R')

tranOffsets = read.csv(allPhenoTranFile,header=TRUE)
tranOffsets$meanDOY <- tranOffsets$meanDOY - 181

#Investigating if data 7 days after might be a better predictor of the inflection point than 7 days before
sitePvals <- matrix(nrow=nrow(siteData),ncol=3)
allRFDat <- matrix(nrow=0,ncol=8)
for(s in (seq_along(siteData$siteName)[-6])){
  siteName <- siteData$siteName[s]
  print(siteName)
  load(paste0(dataDirectory,siteName,"_dataFinal.RData")) #dataFinal <- rescaled p so I also need the original
  load(file=paste0(dataDirectory,siteName,"_phenopixOutputs.RData"))
  
  tranID <- which(tranOffsets$siteName==siteName)
  tran_DOY_mean <- as.numeric(tranOffsets[tranID,2])
  tranDays <- as.numeric(tranOffsets[tranID,3:(2+dataFinal$N)]) + tran_DOY_mean
  
  if(length(dataFinal$years)>1){
    siteSpecifDat <- matrix(nrow=0,ncol=5)
    for(t in seq_along(tranDays)){
      if(!is.na(tranDays[t])){
        for(w in 5:0){
          temp_before <- dataFinal$TairMu[(round(tranDays[t],digits=0)-(7*(w+1))):(round(tranDays[t],digits=0)-(w*7+1))]
          temp_after <- dataFinal$TairMu[(round(tranDays[t],digits=0)-(w*7)):(round(tranDays[t],digits=0)-(w*7)+6)]
          
          d_before <- dataFinal$D[(round(tranDays[t],digits=0)-(7*(w+1))):(round(tranDays[t],digits=0)-(w*7+1))]
          d_after <- dataFinal$D[(round(tranDays[t],digits=0)-(w*7)):(round(tranDays[t],digits=0)-(w*7)+6)]
          siteSpecifDat <- rbind(siteSpecifDat,c(w,mean(temp_before),mean(temp_after),mean(d_before),mean(d_after)))
        }
      }
    }
    siteSpecifDat <- as.data.frame(siteSpecifDat)
    colnames(siteSpecifDat) <- c("w","temp_before","temp_after","d_before","d_after")
    siteSpecifDat$occurred <- (siteSpecifDat$w==0)
    siteSpecifDat$TxD_before <- siteSpecifDat$d_before * siteSpecifDat$temp_before
    siteSpecifDat$TxD_after <- siteSpecifDat$d_after * siteSpecifDat$temp_after
    siteSpecifDat$w <- NULL
    siteSpecifDat$siteName <- siteName
    
    allRFDat <- rbind(allRFDat,siteSpecifDat)
  }
  
}


sumData <- read.csv('allPhenocamDBsitesSummary.csv') %>% dplyr::select(siteName,Lat,MAT)
allRFDat <- left_join(allRFDat,sumData,by="siteName")

j=which(names(allRFDat)=="occurred")
allRFDat$occurred <- as.factor(as.numeric(allRFDat$occurred))
allRFDat$siteName <- as.factor(allRFDat$siteName)
sample = sample.split(allRFDat[,j], SplitRatio = .75)
train = subset(allRFDat, sample == TRUE)
test  = subset(allRFDat, sample == FALSE)

train$siteName <- NULL
rf <- randomForest(occurred ~ ., data=train)

names(rf$importance) <- c("T Before","T After","D Before","D After","T x D Before","T x D After","Lat","MAT")

imp <- importance(rf)#, class=class, scale=scale)#, type=type) 
values <- numeric()
for(i in seq_along(names(imp))){
  values <- c(values,imp[i]) 
} 

jpeg("RF_importancePlot.jpeg",width=8,height=5,units="in",res=1000)
par(las=1,mar=c(5,8,1,1))
barPlot <- barplot(values[order(values)], names.arg=names(imp)[order(values)],horiz = TRUE,
                   xlab="Importance Ranking (Mean Decrease in Gini)",space=20,col="black",xlim=c(0,75))
points(y=barPlot,x=values[order(values)],pch=20,cex=3)
dev.off()