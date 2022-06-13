library('rjags')
library('runjags')
library('ecoforecastR')
library('RColorBrewer')

# sites <- c("harvard","umichbiological","bostoncommon","coweeta","howland2","missouriozarks","queens","dukehw",
#            "lacclair","bbc1","bartlettir","oakridge1","hubbardbrook","readingma","bullshoals",
#            "willowcreek","downerwoods","laurentides","russellsage","sanford","boundarywaters")

# sites <- c("harvard","umichbiological","bostoncommon","coweeta","howland2",
#            "morganmonroe","missouriozarks","queens","dukehw","lacclair","bbc1","NEON.D08.DELA.DP1.00033",
#            "bartlettir","oakridge1","hubbardbrook","alligatorriver","readingma","bullshoals",
#            "willowcreek","downerwoods","laurentides","russellsage","sanford","boundarywaters")

sites <- c("coweeta","howland2","queens","lacclair",
           "bartlettir","hubbardbrook","alligatorriver","readingma",
           "willowcreek","downerwoods","laurentides","sanford","boundarywaters") #1

# sites <- c("umichbiological","morganmonroe","missouriozarks","bbc1",
#            "NEON.D08.DELA.DP1.00033","russellsage") #2
#sites <- "coweeta"

# sites <- c("harvard","coweeta","lacclair",
#            "bartlettir","hubbardbrook","alligatorriver","laurentides")
#sites <- c("howland2", "morganmonroe","alligatorriver","laurentides")

dataDirectory <- "data/"
Nmc <- 1000
#ns <- c(49,56,63,70,77,84,91,98,105,112,119,126,183)
ns <- seq(40,99)
#ns <- seq(35,100)
#ns <- seq(151,183)
#ns <- seq(101,150)
#ns <- seq(35,183)

cols <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f",'#ff7f00','#cab2d6')
NT <- 184
days <- seq(1,NT)
pdfName <- "hindcastsDifferentWindows_forecasts_sub1sub.pdf"

print(pdfName)
pdf(file=pdfName,height=18,width=60)
par(mfrow=c(7,2))

for(s in 1:length(sites)){
  siteName <- sites[s]
  
  print(siteName)
  load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData")) #Load Data
  yearRemoved <- dataFinal$yearRemoved
  
  yearInt <- which(dataFinal$years==yearRemoved)
  for(n in ns){
    fileName <- paste0('finalVarBurns/',siteName,"_meanTemp_summer",n,"_expBreak_slope_forecast_b3_calibration_varBurn.RData")
    if(file.exists(fileName)){
      print(fileName)
      res <- try(load(fileName))#Load Model Output 
      if(inherits(res,"try-error")){
        next
      }
      if(typeof(out.burn)!=typeof(FALSE)){
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
        ecoforecastR::ciEnvelope(seq(1,length(dataFinal$p)),ci[1,],ci[3,],col=col.alpha("blue",0.5))
        points(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,col="red")
        includeP <- dataFinal$p
        includeP[(n+1):nrow(includeP),] <- NA
        includeP[,yearInt] <- NA
        points(seq(1,length(includeP)),includeP,pch=20,col="black")
      }else{
        fileName <- paste0('finalVarBurns/partial_',siteName,"_meanTemp_summer",n,"_expBreak_slope_forecast_b3_calibration_varBurn.RData")
        print(fileName)
        res <- try(load(fileName))#Load Model Output 
        if(inherits(res,"try-error")){
          next
        }
        out.mat <- data.frame(as.matrix(partialOutput$params))
        prow = sample.int(nrow(out.mat),Nmc,replace=TRUE)
        
        b0 <- out.mat$b0[prow]
        b4 <- out.mat$b4[prow]
        b3 <- out.mat$b3[prow]
        pred.mat <- data.frame(as.matrix(partialOutput$predict))
        ci = apply(pred.mat,2,quantile,c(0.025,0.5,0.975))
        
        if(diff(quantile(b0,c(0.025,0.975)))<0.25 & diff(quantile(b4,c(0.025,0.975)))<0.25 & diff(quantile(b3,c(0.025,0.975)))<0.25){
          convergedWell <- TRUE
        }else{
          convergedWell <- FALSE
        }
        
        plot(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,main=paste0(siteName," ",n,convergedWell),ylim=c(0,1),xlim=c(0,2500))
        ecoforecastR::ciEnvelope(seq(1,length(dataFinal$p)),ci[1,],ci[3,],col=col.alpha("green",0.5))
        points(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,col="red")
        includeP <- dataFinal$p
        includeP[(n+1):nrow(includeP),] <- NA
        includeP[,yearInt] <- NA
        points(seq(1,length(includeP)),includeP,pch=20,col="black")
      }
    }else{
      fileName <- paste0('finalVarBurns/partial_',siteName,"_meanTemp_summer",n,"_expBreak_slope_forecast_b3_calibration_varBurn.RData")
      if(file.exists(fileName)){
        print(fileName)
        res <- try(load(fileName))#Load Model Output 
        if(inherits(res,"try-error")){
          next
          }
        out.mat <- data.frame(as.matrix(partialOutput$param))
        prow = sample.int(nrow(out.mat),Nmc,replace=TRUE)
        
        b0 <- out.mat$b0[prow]
        b4 <- out.mat$b4[prow]
        b3 <- out.mat$b3[prow]
        pred.mat <- data.frame(as.matrix(partialOutput$predict))
        ci = apply(pred.mat,2,quantile,c(0.025,0.5,0.975))
        
        if(diff(quantile(b0,c(0.025,0.975)))<0.25 & diff(quantile(b4,c(0.025,0.975)))<0.25 & diff(quantile(b3,c(0.025,0.975)))<0.25){
          convergedWell <- TRUE
        }else{
          convergedWell <- FALSE
        }
        
        plot(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,main=paste0(siteName," ",n,convergedWell),ylim=c(0,1),xlim=c(0,2500))
        ecoforecastR::ciEnvelope(seq(1,length(dataFinal$p)),ci[1,],ci[3,],col=col.alpha("green",0.5))
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
}

dev.off()

