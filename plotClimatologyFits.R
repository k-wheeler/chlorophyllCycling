library('rjags')
library('runjags')
library('ecoforecastR')
library('RColorBrewer')
library(scoringRules)
sites <- c("harvard","umichbiological","bostoncommon","coweeta","howland2",
           "morganmonroe","missouriozarks","queens","dukehw","lacclair","bbc1","NEON.D08.DELA.DP1.00033",
           "arbutuslake","bartlettir","proctor","oakridge1","hubbardbrook","canadaOA","alligatorriver","readingma",
           "bullshoals","thompsonfarm2N","ashburnham","shalehillsczo")
yearsRemoved <- c(2015,2010,2020,2015,2017,2015,2015,2010,2017,2019,2018,2017,
                  2012,2019,2019,2010,2014,2015,2017,2018,2016,2011,2012,2019)
siteYearRemoves <- cbind(sites,yearsRemoved)

dataDirectory <- "data/"
Nmc <- 5000
ns <- c(49,56,63,70,77,84,91,98,105,112,119,126,183)
cols <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f",'#ff7f00','#cab2d6')
NT <- 184
days <- seq(1,NT)
pdfName <- "hindcastsDifferentWindows_forecasts_climatology.pdf"
pdf(file=pdfName,height=18,width=60)
#par(mfrow=c(7,2))
par(mfcol=c(6,2))

for(s in 1:nrow(siteYearRemoves)){
  #for(s in 1:2){
  
  siteName <- siteYearRemoves[s,1]
  
  print(siteName)
  load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData")) #Load Data
  yearRemoved <- as.numeric(siteYearRemoves[siteYearRemoves[,1]==siteName,2])
  
  yearInt <- which(dataFinal$years==yearRemoved)
  climFileName <- paste0(siteName,"_climatology_forecast_calibration_varBurn.RData")
  load(climFileName)
  pred.matYr <- data.frame(as.matrix(out.burn))[,1:184]
  pred.mat <- pred.matYr
  for(yr in 2:dataFinal$N){
    pred.mat <- cbind(pred.mat,pred.matYr)
  }
  ci = apply(pred.mat,2,quantile,c(0.025,0.5,0.975))
  
  plot(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,main=paste0(siteName," ",n),ylim=c(0,1),xlim=c(0,2500))
  ecoforecastR::ciEnvelope(seq(1,length(dataFinal$p)),ci[1,],ci[3,],col=col.alpha("blue",0.5))
  points(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,col="red")
  includeP <- dataFinal$p
  includeP[(n+1):nrow(includeP),] <- NA
  includeP[,yearInt] <- NA
  points(seq(1,length(includeP)),includeP,pch=20,col="black")
}
dev.off()

