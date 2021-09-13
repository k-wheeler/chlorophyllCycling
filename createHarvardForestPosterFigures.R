##Harvard Forest Mid-term Review Poster Figures
library('rjags')
library('runjags')
library('ecoforecastR')
library('RColorBrewer')
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
      mu[i] <- bd[i] + max(syn[i],0)
      mu[i] <- max(0,min(mu[i],IC[min(i,length(IC))]))
      
      if(length(Q)>1){
        xNew <- c(xNew,rnorm(1,mu[i],Q[i]))
      }else{
        xNew <- c(xNew,rnorm(1,mu[i],Q))
      }
    }
    x[,t] <- xNew
    Xprev <- x[,t]
  }
  return(x)
}

load('harvard_b3_final_calibration_varBurn.RData')

Nmc <- 1000
out.mat <- data.frame(as.matrix(out.burn$param))
baseTempOrig <- 20
b0 <- out.mat$b0
b1 <- out.mat$b1
b2 <- out.mat$b2
b3 <- out.mat$b3
prow = sample.int(nrow(out.mat),Nmc,replace=TRUE)

sites <- c("harvard","morganmonroe","asuhighlands","bullshoals")
years <- c(2,1,2,2)
dtes <- seq(as.Date("2021-08-01"),as.Date("2021-12-31"),"day")

pdf(file="UncertaintyAnalysis_HFposterFigures.pdf",height=6,width=8)
for(s in 1:length(sites)){
  fileName <- paste0(sites[s],"_b3_final_calibration_varBurn.RData")
  print(fileName)
  tle <- strsplit(fileName,"_")[[1]][1:2]
  siteName <- strsplit(fileName,"_")[[1]][1]
  sites <- c(sites,siteName)
  load(fileName)
  load(paste(siteName,"_CDD",baseTempOrig,"_dataFinal.RData",sep=""))
  out.mat <- data.frame(as.matrix(out.burn$param))
  
  NT <- length(dataFinal$TairMu[,years[s]])
  out.mat.pred <- data.frame(as.matrix(out.burn$predict))
  
  prow.pred = sample.int(nrow(out.mat),Nmc,replace=TRUE)
  initialXs <- out.mat.pred[prow.pred,years[s]]
  days <- seq(1,NT)
  ysDet <- forecastStep(IC=mean(initialXs),b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),n=1,NT=length(days),Tair=dataFinal$TairMu[,years[s]]) #Dependent on IC, but independent for >=1
  ysParam <- forecastStep(IC=mean(initialXs),b0=b0[prow],b1=b1[prow],b2=b2[prow],b3=b3[prow],n=Nmc,NT=NT,Tair=dataFinal$TairMu[,years[s]])
  
  ysProc <- forecastStep(IC=mean(initialXs),b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),n=Nmc,NT=NT,Q=sqrt(1/out.mat$p.proc[prow.pred]),Tair=dataFinal$TairMu[,years[s]])
  ysIC <- forecastStep(IC=initialXs,b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),n=length(initialXs),NT=NT,Tair=dataFinal$TairMu[,years[s]])
  
  N.IP.ci = apply(ysParam,2,quantile,c(0.025,0.5,0.975))
  N.Proc.ci = apply(ysProc,2,quantile,c(0.025,0.5,0.975))
  N.IC.ci = apply(ysIC,2,quantile,c(0.025,0.5,0.975))

  plot(dtes,dataFinal$p[,years[s]],pch=20,main=paste(siteName,years[s]),ylim=c(0,1),xlab="Time",ylab="Rescaled GCC")
  
  ecoforecastR::ciEnvelope(dtes,N.Proc.ci[1,],N.Proc.ci[3,],col=col.alpha("#66c2a5",1)) #green
  ecoforecastR::ciEnvelope(dtes,N.IP.ci[1,],N.IP.ci[3,],col=col.alpha("#fc8d62",1)) #orange
  ecoforecastR::ciEnvelope(dtes,N.IC.ci[1,],N.IC.ci[3,],col=col.alpha("#8da0cb",1)) #purple 
  points(dtes,dataFinal$p[,years[s]],pch=20)
  lines(dtes,ysDet,col="black",lwd=3)

}
plot(numeric(),numeric(),ylim=c(0,1),xlim=c(0,1))
legend("topleft",col=c("black","black","#8da0cb","#fc8d62","#66c2a5"),
       c("PhenoCam Observations","Mean Estimate","Initial Condition Error","Parameter Error",
         "Process Error"),pch=c(20,15,15,15,15))
dev.off()

