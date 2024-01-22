#Investigate greenness, temperature, and photoperiod around transitions 
library(scales)
library('mgcv')
library(tidyverse)
library(randomForest)
library(caTools)
library(tidyverse)
dataDirectory <- "data/"
siteData <- read.csv('/projectnb/dietzelab/kiwheel/chlorophyllCycling/allPhenocamDBsitesComplete.csv',header=TRUE)

tranOffsets = read.csv('phenocamTransitions_fromMean.csv',header=TRUE)
tranOffsets$meanDOY <- tranOffsets$meanDOY - 181

greenData <- matrix(nrow=nrow(siteData),ncol=7)
coData <- matrix(nrow=nrow(siteData),ncol=7)
#slopeData <- matrix(nrow=nrow(siteData),ncol=3)

for(s in seq_along(siteData$siteName)[-6]){
  siteName <- siteData$siteName[s]
  print(siteName)
  load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData")) #dataFinal <- rescaled p so I also need the original
  load(file=paste0(dataDirectory,siteName,"_phenopixOutputs.RData"))
  
  tranID <- which(tranOffsets$siteName==siteName)
  tran_DOY_mean <- as.numeric(tranOffsets[tranID,2])
  tranDays <- as.numeric(tranOffsets[tranID,3:(2+dataFinal$N)]) + tran_DOY_mean
  
  if(length(dataFinal$years)>1){
    
    rescalingMin <- mean(allDat[,2],na.rm=TRUE)
    rescalingMax <- mean(allDat[,3],na.rm=TRUE)
    dataFinal$gcc <- rescale(dataFinal$p,to=c(rescalingMin,rescalingMax)) #Not quite exact
    p_trans <- numeric()
    gcc_trans <- numeric()
    temp_before <- numeric()
    temp_after <- numeric()
    temp_diff <- numeric()
    for(t in seq_along(tranDays)){
      p_trans <- c(p_trans,dataFinal$p[round(tranDays[t],digits=0),t])
      gcc_trans <- c(gcc_trans,dataFinal$gcc[round(tranDays[t],digits=0),t])
      # temp_before <- c(temp_before,mean(dataFinal$TairMu[((round(tranDays[t],digits=0)-7)):(round(tranDays[t],digits=0)-1),t]))
      # temp_after <- c(temp_after,mean(dataFinal$TairMu[((round(tranDays[t],digits=0))):(round(tranDays[t],digits=0)+6),t]))
      # temp_diff <- c(temp_diff,(mean(dataFinal$TairMu[((round(tranDays[t],digits=0))):(round(tranDays[t],digits=0)+6),t])-
      #                             mean(dataFinal$TairMu[((round(tranDays[t],digits=0)-7)):(round(tranDays[t],digits=0)-1),t])))
      
      temp_before <- c(temp_before,mean(dataFinal$TairMu[((round(tranDays[t],digits=0)-14)):(round(tranDays[t],digits=0)-8),t]))
      temp_after <- c(temp_after,mean(dataFinal$TairMu[((round(tranDays[t],digits=0)-7)):(round(tranDays[t],digits=0)-1),t]))
      temp_diff <- c(temp_diff,(mean(dataFinal$TairMu[((round(tranDays[t],digits=0)-7)):(round(tranDays[t],digits=0)-1),t])-
                                  mean(dataFinal$TairMu[((round(tranDays[t],digits=0)-14)):(round(tranDays[t],digits=0)-8),t])))
      beforeGCC <- dataFinal$gcc[1:(round(tranDays[t],digits=0)-1),t]
      mdl_g <- lm(beforeGCC ~ seq_along(beforeGCC))
      print(mdl_g$coefficients[2])
      # beforeP <- dataFinal$p[1:(round(tranDays[t],digits=0)-1),t]
      # mdl_p <- lm(beforeP ~ seq_along(beforeP))
      
      # temp_before <- c(temp_before,mean(dataFinal$TairMu[((round(tranDays[t],digits=0)-14)):(round(tranDays[t],digits=0)-8),t]))
      # temp_after <- c(temp_after,mean(dataFinal$TairMu[((round(tranDays[t],digits=0))-7):(round(tranDays[t],digits=0)-1),t]))
    }
    
    greenData[s,] <- c(siteName,dataFinal$N,mean(p_trans,na.rm=TRUE),sd(p_trans,na.rm=TRUE),mean(gcc_trans,na.rm=TRUE),sd(gcc_trans,na.rm=TRUE),
                       rescalingMax-rescalingMin)
    coData[s,] <- c(siteName,dataFinal$N,mean(temp_before),sd(temp_before),mean(temp_after),sd(temp_after),(sum(temp_diff<0)/dataFinal$N))
  }
}

greenData <- as.data.frame(greenData)
colnames(greenData) <- c('siteName','numYears',"p_mean","p_sd","gcc_mean","gcc_sd","gcc_range")  
greenData$gccSDPercRange <- as.numeric(greenData$gcc_sd)/as.numeric(greenData$gcc_range)*100

coData <- as.data.frame(coData)
colnames(coData) <- c('siteName','numYears',"tempBefore_mean","tempBefore_sd","tempAfter_mean","tempAfter_sd","percentYearsDecrease") 
coData$tempDiff <- as.numeric(coData$tempAfter_mean)-as.numeric(coData$tempBefore_mean)
sum(na.omit(coData$tempDiff)<0)/(nrow(coData)-1)

# plot(density(greenData$gccSDPercRange))
# plot(density(as.numeric(na.omit(coData$percentYearsDecrease))))
# lines(density(as.numeric(na.omit(coData$percentYearsDecrease))),col="green")

#Investigating if data 7 days after might be a better predictor of the inflection point than 7 days before
sitePvals <- matrix(nrow=nrow(siteData),ncol=3)
allRFDat <- matrix(nrow=0,ncol=8)
for(s in (seq_along(siteData$siteName)[-6])){
  siteName <- siteData$siteName[s]
  print(siteName)
  load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData")) #dataFinal <- rescaled p so I also need the original
  load(file=paste0(dataDirectory,siteName,"_phenopixOutputs.RData"))
  
  tranID <- which(tranOffsets$siteName==siteName)
  tran_DOY_mean <- as.numeric(tranOffsets[tranID,2])
  tranDays <- as.numeric(tranOffsets[tranID,3:(2+dataFinal$N)]) + tran_DOY_mean
  
  if(length(dataFinal$years)>1){
    siteSpecifDat <- matrix(nrow=0,ncol=5)
    for(t in seq_along(tranDays)){
      for(w in 5:0){
        temp_before <- dataFinal$TairMu[(round(tranDays[t],digits=0)-(7*(w+1))):(round(tranDays[t],digits=0)-(w*7+1))]
        temp_after <- dataFinal$TairMu[(round(tranDays[t],digits=0)-(w*7)):(round(tranDays[t],digits=0)-(w*7)+6)]
        
        d_before <- dataFinal$D[(round(tranDays[t],digits=0)-(7*(w+1))):(round(tranDays[t],digits=0)-(w*7+1))]
        d_after <- dataFinal$D[(round(tranDays[t],digits=0)-(w*7)):(round(tranDays[t],digits=0)-(w*7)+6)]
        #prod_before <- temp_before * d_before
        #prod_after <- temp_after * temp_after
        siteSpecifDat <- rbind(siteSpecifDat,c(w,mean(temp_before),mean(temp_after),mean(d_before),mean(d_after)))
      }
    }
    siteSpecifDat <- as.data.frame(siteSpecifDat)
    colnames(siteSpecifDat) <- c("w","temp_before","temp_after","d_before","d_after")
    siteSpecifDat$occurred <- (siteSpecifDat$w==0)
    siteSpecifDat$TxD_before <- siteSpecifDat$d_before * siteSpecifDat$temp_before
    siteSpecifDat$TxD_after <- siteSpecifDat$d_after * siteSpecifDat$temp_after
    # yesData <- siteSpecifDat %>% filter(occurred)
    # noData <- siteSpecifDat %>% filter(!occurred)
    
    siteSpecifDat$w <- NULL
    siteSpecifDat$siteName <- siteName
    # j=which(names(siteSpecifDat)=="occurred")
    # sample = sample.split(siteSpecifDat[,j], SplitRatio = .75)
    # train = subset(siteSpecifDat, sample == TRUE)
    # 
    # rf <- randomForest(occurred ~ ., data=train)
    # 
    # sitePvals[s,] <- c(siteName,t.test(yesData$prod_before,noData$prod_before)$p.value,
    #                    t.test((yesData$prod_after/yesData$prod_before),(noData$prod_after/noData$prod_before))$p.value))
    # 
    # 
    # mdl_before <- gam(occurred ~ s(temp_before,k=5)+s(d_before,k=5),family = binomial(link="logit"),
    #                   data=siteSpecifDat,
    #                   method = "REML")
    # mdl_before <- lm(factor(occurred)~temp_before +d_before + temp_before * d_before,
    #                  data=siteSpecifDat)
    # sm <- summary(mdl_after)
    # r_before <- sm$r.sq
    # 
    # mdl_after <- gam(occurred ~ s(temp_after,k=5)+s(d_after,k=5),family = binomial(link="logit"),
    #                   data=siteSpecifDat,
    #                   method = "REML")
    # sm <- summary(mdl_after)
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

#partialPlot(x=rf,pred.data=train,x.var="prod_after")

yesData <- allRFDat %>% filter(occurred==1)
#plot(density(as.numeric(yesData$TxD_after)))

names(rf$importance) <- c("T Before","T After","D Before","D After","T x D Before","T x D After","Lat","MAT")

imp <- importance(rf, class=class, scale=scale, type=type)
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

#randomForest::varImpPlot(rf,main = "Random Forest Importance Ranking") #KWheeler Edited

test_x <- data.matrix(test[,-j])
test_y <- data.matrix(test[,j])
pred_y <- predict(rf,test_x)
sum(pred_y==test_y)/length(pred_y) #Percent that we predicted correctly

sitePvals <- as.data.frame(sitePvals)
sitePvals$Decrease<- sitePvals[,3]<sitePvals[,2]
sum(na.omit(sitePvals$Decrease))/length(na.omit(sitePvals$Decrease))
