library('rjags')
library('runjags')
library('ecoforecastR')
library('RColorBrewer')

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
pdfName <- "hindcastsDifferentWindows_forecasts.pdf"
pdf(file=pdfName,height=18,width=60)
par(mfrow=c(7,2))

for(s in 1:nrow(siteYearRemoves)){
  allb0s <- matrix(ncol=0,nrow=Nmc)
  allb3s <- matrix(ncol=0,nrow=Nmc)
  allb4s <- matrix(ncol=0,nrow=Nmc)
  
  siteName <- siteYearRemoves[s,1]
  
  print(siteName)
  load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData")) #Load Data
  yearRemoved <- as.numeric(siteYearRemoves[siteYearRemoves[,1]==siteName,2])
  
  yearInt <- which(dataFinal$years==yearRemoved)
  for(n in ns){
    fileName <- paste0(siteName,"_meanTemp_summer",n,"_expBreak_slope_forecast_b3_calibration_varBurn.RData")
    if(file.exists(fileName)){
      print(fileName)
      load(fileName) #Load Model Output 
      out.mat <- data.frame(as.matrix(out.burn$param))
      prow = sample.int(nrow(out.mat),Nmc,replace=TRUE)
      
      b0 <- out.mat$b0[prow]
      b4 <- out.mat$b4[prow]
      b3 <- out.mat$b3[prow]
      pred.mat <- data.frame(as.matrix(out.burn$predict))
      ci = apply(pred.mat,2,quantile,c(0.025,0.5,0.975))
      
      plot(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,main=paste0(siteName," ",n),ylim=c(0,1),xlim=c(0,2500))
      ecoforecastR::ciEnvelope(seq(1,length(dataFinal$p)),ci[1,],ci[3,],col=col.alpha("blue",0.5))
      points(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,col="red")
      includeP <- dataFinal$p
      includeP[(n+1):nrow(includeP),] <- NA
      includeP[,yearInt] <- NA
      points(seq(1,length(includeP)),includeP,pch=20,col="black")
      
    }else{
      fileName <- paste0('partial_',siteName,"_meanTemp_summer",n,"_expBreak_slope_forecast_b3_calibration_varBurn.RData")
      if(file.exists(fileName)){
        print(fileName)
        load(fileName) #Load Model Output 
        out.mat <- data.frame(as.matrix(partialOutput$param))
        prow = sample.int(nrow(out.mat),Nmc,replace=TRUE)
        
        b0 <- out.mat$b0[prow]
        b4 <- out.mat$b4[prow]
        b3 <- out.mat$b3[prow]
        pred.mat <- data.frame(as.matrix(partialOutput$predict))
        ci = apply(pred.mat,2,quantile,c(0.025,0.5,0.975))
        
        plot(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,main=paste0(siteName," ",n),ylim=c(0,1),xlim=c(0,2500))
        ecoforecastR::ciEnvelope(seq(1,length(dataFinal$p)),ci[1,],ci[3,],col=col.alpha("green",0.5))
        points(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,col="red")
        includeP <- dataFinal$p
        includeP[(n+1):nrow(includeP),] <- NA
        includeP[,yearInt] <- NA
        points(seq(1,length(includeP)),includeP,pch=20,col="black")
      }else{
        plot(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,main=paste0(siteName," ",n," (Removed: ",yearInt,")"),ylim=c(0,1),xlim=c(0,2500))
        b0 <- rep(NA,Nmc)
        b3 <- rep(NA,Nmc)
        b4 <- rep(NA,Nmc)
      }
    }
    allb0s <- cbind(allb0s,b0)
    allb3s <- cbind(allb3s,b3)
    allb4s <- cbind(allb4s,b4)
  }
  plot(numeric(),numeric(),ylim=c(0,100),xlim=c(-1,0),main="b0")
  for(c in 1:ncol(allb0s)){
    if(!is.na(allb0s[1,c])){
      lines(density(allb0s[,c]),col=cols[c])
    }
  }
  # plot(numeric(),numeric(),ylim=c(0,10000),xlim=c(0,0.1),main="b3")
  # for(c in 1:ncol(allb3s)){
  #   if(!is.na(allb3s[1,c])){
  #     lines(density(allb3s[,c]),col=cols[c])
  #   }
  # }
  # plot(numeric(),numeric(),ylim=c(0,100),xlim=c(-1,0),main="b4")
  # for(c in 1:ncol(allb4s)){
  #   if(!is.na(allb4s[1,c])){
  #     lines(density(allb4s[,c]),col=cols[c])
  #   }
  # }
  legend('topleft',lty=rep(1,length(cols)),col=cols,legend=ns)
}
dev.off()

