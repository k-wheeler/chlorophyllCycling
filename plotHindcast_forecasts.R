library('rjags')
library('runjags')
library('RColorBrewer')
source('generalVariables.R')
source('ciEnvelope.R')

Nmc <- 1000
ns <- seq(40,183)

cols <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f",'#ff7f00','#cab2d6')
NT <- 184
days <- seq(1,NT)
pdfName <- "hindcastsDifferentWindows_forecasts.pdf"

print(pdfName)
pdf(file=pdfName,height=18,width=60)
par(mfrow=c(7,2))

for(s in 1:length(sites)){
  siteName <- sites[s]
  
  print(siteName)
  load(paste0(dataDirectory,siteName,"_dataFinal.RData")) #Load Data
  yearRemoved <- dataFinal$yearRemoved
  
  yearInt <- which(dataFinal$years==yearRemoved)
  for(n in ns){
    fileName <- paste0(CCmodelOutputsFolder,siteName,"_",n,"_ccModel_forecast_calibration_varBurn.RData")
    if(file.exists(fileName)){
      print(fileName)
      res <- try(load(fileName))#Load Model Output 
      if(inherits(res,"try-error")){
        next
      }
      out.mat <- data.frame(as.matrix(out.burn$params))
      prow = sample.int(nrow(out.mat),Nmc,replace=TRUE)
      
      b0 <- out.mat$b0[prow]
      b4 <- out.mat$b4[prow]
      b3 <- out.mat$b3[prow]
      pred.mat <- data.frame(as.matrix(out.burn$predict))
      ci = apply(pred.mat,2,quantile,c(0.025,0.5,0.975))
      
      if(diff(quantile(b0,c(0.025,0.975)))<0.25 & diff(quantile(b4,c(0.025,0.975)))<0.25 & diff(quantile(b3,c(0.025,0.975)))<0.25){
        convergedWell <- TRUE
      }else{
        convergedWell <- FALSE
      }
      
      plot(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,main=paste0(siteName," ",n,convergedWell),ylim=c(0,1),xlim=c(0,2500))
      ciEnvelope(seq(1,length(dataFinal$p)),ci[1,],ci[3,],col=col.alpha("blue",0.5))
      points(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,col="red")
      includeP <- dataFinal$p
      includeP[(n+1):nrow(includeP),] <- NA
      includeP[,yearInt] <- NA
      points(seq(1,length(includeP)),includeP,pch=20,col="black")
    }else{
      plot(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,main=paste0(siteName," ",n," (Removed: ",yearInt,")"),ylim=c(0,1),xlim=c(0,2500))
    }
  }
}

dev.off()

