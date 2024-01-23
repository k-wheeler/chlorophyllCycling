###Paper Figures
library('rjags')
library('runjags')
library('RColorBrewer')
library('scoringRules')
source('generalVariables')

tranOffsets <- read.csv(allPhenoTranFile,header=TRUE)
NT <- 184
days <- seq(1,NT)
crpsDiffsALL <- matrix(ncol=187,nrow=0)
daysALL <- matrix(ncol=187,nrow=0)
allTrans <- numeric()
Nmc=1000
for(s in 1:length(sites)){
  siteName <- sites[s]
  print(siteName)
  load(paste0(dataDirectory,siteName,"_dataFinal.RData")) #Load Data
  yearRemoved <- dataFinal$yearRemoved

  yearInt <- which(dataFinal$years==yearRemoved)
  climFileName <- paste0(climatologyModelOutputsFolder,siteName,"_climatology_forecast_calibration_varBurn.RData")
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
  for(n in 35:183){
    print(n)
    fileName <- paste0(CCmodelOutputsFolder,siteName,"_",n,"_ccModel_forecast_calibration_varBurn.RData")
    if(file.exists(fileName)){
      res <- try(load(fileName))#Load Model Output
      if(inherits(res,"try-error")){
        next
      }
    }else{
      next #No files
    }

    out.mat <- data.frame(as.matrix(out.burn$param))
    b0 <- out.mat$b0
    b4 <- out.mat$b4
    b3 <- out.mat$b3
    pred.mat <- data.frame(as.matrix(out.burn$predict))
    prow = sample.int(nrow(pred.mat),Nmc,replace=TRUE)

    if(diff(quantile(b0,c(0.025,0.975)))<0.25 & diff(quantile(b4,c(0.025,0.975)))<0.25 & diff(quantile(b3,c(0.025,0.975)))<0.25){
      convergedWell <- TRUE
    }else{
      convergedWell <- FALSE
    }

    for(yr in 1:length(dataFinal$years)){
      yrName <- dataFinal$years[yr]
      pred.mat.year <- pred.mat[prow,(184*(yr-1)+1):(184*(yr-1)+184)]

      crpsVals <- numeric()
      for(i in 1:ncol(pred.mat.year)){
        if(!is.na(dataFinal$p[i,yr])){
          crpsVals <- c(crpsVals,crps_sample(y=dataFinal$p[i,yr],dat=pred.mat.year[,i]))
        }else{
          crpsVals <- c(crpsVals,NA)
        }
      }
      crps_diff <- crpsVals-(crpsClim[(184*(yr-1)+1):(184*(yr-1)+184)])
      tranID <- which(tranOffsets$siteName==siteName)
      tran_DOY <- as.numeric(tranOffsets[tranID,(yr+2)]) + as.numeric(tranOffsets[tranID,2])

      offset <- round(((n+182) - tran_DOY),digits=0)
      if(convergedWell){
        crpsDiffsALL <- rbind(crpsDiffsALL,c(siteName,yrName,offset,crps_diff))
        daysALL <- rbind(daysALL,c(siteName,yrName,NA,((182 +days)-tran_DOY)))
        allTrans <- c(allTrans,as.numeric(tranOffsets[tranID,(yr+2)]))
      }
    }
  }
}

daysOffset <- seq(round(min((as.numeric(daysALL[,4:ncol(daysALL)]))),digits=0),
                  round(max((as.numeric(daysALL[,4:ncol(daysALL)]))),digits=0))
reorganizedDat <- data.frame(matrix(nrow=nrow(daysALL),ncol=3+length(daysOffset)))
reorganizedDat[,1] <- crpsDiffsALL[,1]
reorganizedDat[,2] <- crpsDiffsALL[,2]
reorganizedDat[,3] <- crpsDiffsALL[,3]


for(i in 1:nrow(reorganizedDat)){ #rows mean different forecasts
  for(d in 1:length(daysOffset)){
    ind <- which(round(as.numeric(daysALL[i,]),digits=0)==as.numeric(daysOffset[d])) #forecasted day
    if(length(ind)>0){
      ind <- ind[1]
      reorganizedDat[i,(d+3)] <- as.numeric(crpsDiffsALL[i,ind])
    }
  }
}
save(reorganizedDat,file=reorganizedDatFile)
save(daysOffset,file=daysOffsetFile)

##Pad reorganizedDat to replace missing offset values at the end with last available data
for(s in 1:length(sites)){
  siteDat <- reorganizedDat[reorganizedDat[,1]==sites[s],]
  siteYears <- as.numeric(unique(siteDat[,2]))
  for(yr in 1:length(siteYears)){
    siteYearDat <- siteDat[siteDat[,2]==siteYears[yr],]
    for(n in daysOffset){
      if((!(n %in% siteYearDat[,3]))& (n>as.numeric(siteYearDat[nrow(siteYearDat),3]))){
        newRow <- siteYearDat[nrow(siteYearDat),]
        newRow[3] <- n
        reorganizedDat <- rbind(reorganizedDat,newRow)
      }
    }
  }
}

crpsMat <- matrix(ncol=length(daysOffset),nrow=length(daysOffset))

for(cl in 1:length(daysOffset)){ #the columns
  print(cl)
  for(d in 1:length(daysOffset)){ #the values
    subDat <- reorganizedDat[as.numeric(reorganizedDat$X3)==as.numeric(daysOffset[cl]),] #Which forecasts
    if(nrow(subDat)>0){
      crpsMat[d,cl] <- sum(na.omit(subDat[,(d+3)])<0)/length(na.omit(subDat[,(d+3)]))
    }
  }
}
save(crpsMat,file=crpsMatFile)

image(x=daysOffset,y=daysOffset,z=t(crpsMat),col=c("#a6611a","#dfc27d","#80cdc1","#018571"),xlab="Days Included Relative to Transition Date",
      ylab="Predicted Day Relative to Transition Date",main="All")
