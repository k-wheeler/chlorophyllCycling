library('rjags')
library('runjags')
library('ecoforecastR')
library('RColorBrewer')

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

#Create forecast step model:
#' Basic logistic forecast step
#'
#' @param IC Initial conditions
#' @param b0 The parameter b0
#' @param b1 The parameter b1
#' @param b2 The parameter b2
#' @param b3 The parameter b3
#' @param Q Process error (default = 0 for deterministic runs)
#' @param n Size of Monte Carlo ensemble
#' @param NT number of days for forecast
#' @import rjags
#' @import runjags
#' @import ecoforecastR
#'
#' @return
#' @export
#'
#' @examples
forecastStep <- function(IC,b0,b1,b2,b3,Q=0,n,NT,Tair){
  x <- matrix(NA,n,NT)
  if(length(IC)==1){
    Xprev <- rep(IC,n)
  }else{
    Xprev <- IC
  }
  
  for(t in 1:NT){
    bd <- Xprev + b0 + (b1 * Xprev) + (b2 * Xprev **2) 
    syn <- b3 * Tair[t]
    if(length(syn)==1){
      syn <- rep(syn,n)
    }
    
    xNew <- numeric()
    mu <- rep(NA,n)
    for(i in 1:n){
      #mu[i] <- max(0,min(mu[i],0.999))
      mu[i] <- bd[i] + max(syn[i],0)
      mu[i] <- max(0,min(mu[i],IC[min(i,length(IC))]))
      
      if(length(Q)>1){
        xNew <- c(xNew,rnorm(1,mu[i],Q[i]))
      }else{
        xNew <- c(xNew,rnorm(1,mu[i],Q))
      }
      #xNew <- c(xNew,max(0, min(0.99,xl)))
    }
    x[,t] <- xNew ## trunate normal process error
    Xprev <- x[,t]
  }
  return(x)
}

load('harvard_b3_final_calibration_varBurn.RData')

Nmc <- 1000
out.mat <- data.frame(as.matrix(out.burn$param))

b0 <- out.mat$b0
b1 <- out.mat$b1
b2 <- out.mat$b2
b3 <- out.mat$b3
prow = sample.int(nrow(out.mat),Nmc,replace=TRUE)

files <- dir(pattern = "_b3_final_calibration_varBurn.RData")

sites <- character()

pdf(file="UncertaintyAnalysis_noTrigger_b3_calibration_harvardParas.pdf",height=6,width=40)
for(f in 1:length(files)){
  
  fileName <- files[f]
  if(strsplit(fileName,"_")[[1]][1]!="partial"){
    print(fileName)
    tle <- strsplit(fileName,"_")[[1]][1:2]
    siteName <- strsplit(fileName,"_")[[1]][1]
    sites <- c(sites,siteName)
    
    print(siteName)
    if(TRUE){
      
      baseTempOrig <- 20
      load(fileName)
      load(paste(siteName,"_CDD",baseTempOrig,"_dataFinal.RData",sep=""))
      if(ncol(dataFinal$TairMu)>1){
        out.mat <- data.frame(as.matrix(out.burn$param))
        
        NT <- length(dataFinal$TairMu[,1])
        
        out.mat.pred <- data.frame(as.matrix(out.burn$predict))
        
        
        prow.pred = sample.int(nrow(out.mat),Nmc,replace=TRUE)
        initialXs <- out.mat.pred[prow.pred,1]
        print(mean(initialXs))
        
        days <- seq(1,NT)
        ysDet <- forecastStep(IC=mean(initialXs),b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),n=1,NT=length(days),Tair=dataFinal$TairMu[,1]) #Dependent on IC, but independent for >=1
        ysParam <- forecastStep(IC=mean(initialXs),b0=b0[prow],b1=b1[prow],b2=b2[prow],b3=b3[prow],n=Nmc,NT=NT,Tair=dataFinal$TairMu[,1])
        
        ysProc <- forecastStep(IC=mean(initialXs),b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),n=Nmc,NT=NT,Q=sqrt(1/out.mat$p.proc[prow.pred]),Tair=dataFinal$TairMu[,1])
        ysIC <- forecastStep(IC=initialXs,b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),n=length(initialXs),NT=NT,Tair=dataFinal$TairMu[,1])
        
        for(i in 2:dataFinal$N){
          print(colnames(out.mat.pred)[dataFinal$n*(i-1)+1])
          initialXs <- out.mat.pred[prow.pred,dataFinal$n*(i-1)+1]
          
          ysDet <- cbind(ysDet,forecastStep(IC=mean(initialXs),b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),n=1,NT=length(days),Tair=dataFinal$TairMu[,i]))
          ysParam <- cbind(ysParam,forecastStep(IC=mean(initialXs),b0=b0[prow],b1=b1[prow],b2=b2[prow],b3=b3[prow],n=Nmc,NT=NT,Tair=dataFinal$TairMu[,i]))
          ysProc <- cbind(ysProc,forecastStep(IC=mean(initialXs),b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),n=Nmc,NT=NT,Q=sqrt(1/out.mat$p.proc[prow.pred]),Tair=dataFinal$TairMu[,i]))
          ysIC <- cbind(ysIC,forecastStep(IC=initialXs,b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),n=length(initialXs),NT=NT,Tair=dataFinal$TairMu[,i]))
        }
        N.IP.ci = apply(ysParam,2,quantile,c(0.025,0.5,0.975))
        N.Proc.ci = apply(ysProc,2,quantile,c(0.025,0.5,0.975))
        N.IC.ci = apply(ysIC,2,quantile,c(0.025,0.5,0.975))
        
        plot(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,main=paste(tle, "p"),ylim=c(0,1))
        
        ecoforecastR::ciEnvelope(seq(1,length(dataFinal$p)),N.Proc.ci[1,],N.Proc.ci[3,],col=col.alpha("green",0.5))
        ecoforecastR::ciEnvelope(seq(1,length(dataFinal$p)),N.IP.ci[1,],N.IP.ci[3,],col=col.alpha("blue",0.5))
        ecoforecastR::ciEnvelope(seq(1,length(dataFinal$p)),N.IC.ci[1,],N.IC.ci[3,],col=col.alpha("pink",0.5))
        lines(seq(1,length(dataFinal$p)),N.Proc.ci[2,],col="green",lwd=2)
        lines(seq(1,length(dataFinal$p)),N.IP.ci[2,],col="blue",lwd=2)
        lines(seq(1,length(dataFinal$p)),N.IC.ci[2,],col="pink",lwd=2)
        lines(seq(1,length(dataFinal$p)),ysDet,col="cyan",lwd=3)
        
        print("done second plot")
        
      }
    }
  }
}

dev.off()