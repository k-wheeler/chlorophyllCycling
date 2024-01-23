library('rjags')
library('runjags')
library('RColorBrewer')
library('scoringRules')
source('generalVariables.R')
source('ciEnvelope.R')

Nmc <- 5000

NT <- 184
totalSiteYears <- 0

pdfName <- "climatology_fits.pdf"
pdf(file=pdfName,height=18,width=60)

par(mfcol=c(6,2))

for(s in 1:nrow(siteData)){
  siteName <- as.character(siteData$siteName[s])
  
  print(siteName)
  load(paste0(dataDirectory,siteName,"_dataFinal.RData")) #Load Data
  totalSiteYears <- totalSiteYears + dataFinal$N
  
  yearInt <- which(dataFinal$years==dataFinal$yearRemoved)
  climFileName <- paste0(climatologyModelOutputsFolder,siteName,"_climatology_forecast_calibration_varBurn.RData")
  if(file.exists(climFileName)){
    load(climFileName)
    pred.matYr <- data.frame(as.matrix(out.burn))[,1:184]
    pred.mat <- pred.matYr
    for(yr in 2:dataFinal$N){
      pred.mat <- cbind(pred.mat,pred.matYr)
    }
    ci = apply(pred.mat,2,quantile,c(0.025,0.5,0.975))
    
    plot(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,main=paste0(siteName," Climatology"),ylim=c(0,1),xlim=c(0,2500),xlab="Day",ylab="GCC")
    ciEnvelope(seq(1,length(dataFinal$p)),ci[1,],ci[3,],col=col.alpha("blue",0.5))
    points(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,col="red")
    includeP <- dataFinal$p
    includeP[,yearInt] <- NA
    points(seq(1,length(includeP)),includeP,pch=20,col="black")
  }else{
    plot(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,main=paste0(siteName," Climatology"),ylim=c(0,1),xlim=c(0,2500),xlab="Day",ylab="GCC")
    points(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,col="red")
    includeP <- dataFinal$p
    includeP[,yearInt] <- NA
    points(seq(1,length(includeP)),includeP,pch=20,col="black")
  }
}
dev.off()

