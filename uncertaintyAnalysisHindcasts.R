library('rjags')
library('runjags')
library('ecoforecastR')
library('RColorBrewer')
library(doParallel)

n.cores <- 24

#register the cores.
registerDoParallel(cores=n.cores)


#Create forecast step model:
#' Basic logistic forecast step
#'
#' @param IC Initial conditions
#' @param b0 The parameter b0
#' @param b1 The parameter b1
#' @param b2 The parameter b2
#' @param b3 The parameter b3
#' @param b4 The parameter b4
#' @param b5 The parameter b5
#' @param Q Process error (default = 0 for deterministic runs)
#' @param n Size of Monte Carlo ensemble
#' @param NT number of days for forecast
#' @param Tair
#' @param D
#' @import rjags
#' @import runjags
#' @import ecoforecastR
#'
#' @return
#' @export
#'
#' @examples
forecastStep <- function(IC,b0,b1,b2,b3,b4,b5,Q=0,n,NT,Tair,D){
  x <- matrix(NA,n,NT)
  if(length(IC)==1){
    Xprev <- rep(IC,n)
  }else{
    Xprev <- IC
  }
  
  for(t in 1:NT){
    bd <- Xprev + b0 + (b1 * Xprev) + (b2 * Xprev **2) 
    syn <- b3 * Tair[t] + b4 * Tair[t] * D[t] + b5 * D[t]
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

uncertaintyAnalysisHindcast <- function(filePattern,pdfName,b,tID){
  
  dataDirectory <- "data/"
  siteData <- read.csv('/projectnb/dietzelab/kiwheel/chlorophyllCycling/allPhenocamDBsitesComplete.csv',header=TRUE)
  
  files <- dir(pattern = filePattern)
  foreach(s=1:nrow(siteData)) %dopar% {
    siteName <- as.character(siteData$siteName[s])
    
    load(paste0(dataDirectory,siteName,"_dataFinal.RData")) #Load Data
    
    for(f in 1:length(files)){
      fileName <- files[f]
      if(strsplit(fileName,"_")[[1]][1]!="partial"){
        print(fileName)
        tle <- strsplit(fileName,"_")[[1]][1:2]
        calibrationSiteName <- strsplit(fileName,"_")[[1]][1]
        outFileName <- paste0("hindcasts/",siteName,"_",calibrationSiteName,"_parameters_",tID,"_",b,"_hindcasts.RData")
        if(!file.exists(outFileName)){
          
          load(fileName) #Load Model Output 
          
          out.mat <- data.frame(as.matrix(out.burn$param))
          
          b0 <- out.mat$b0
          b1 <- out.mat$b1
          b2 <- out.mat$b2
          if(b=="b3"){
            b3 <- out.mat$b3
            b4 <- rep(0,length(b3))
            b5 <- rep(0,length(b3))
          }else if(b=="b4"){
            b4 <- out.mat$b4
            b3 <- rep(0,length(b4))
            b5 <- rep(0,length(b4))
          }else if(b=="b5"){
            b5 <- out.mat$b5
            b3 <- rep(0,length(b5))
            b4 <- rep(0,length(b5))
          }
          
          NT <- length(dataFinal$TairMu[,1])
          
          Nmc <- 1000
          prow = sample.int(nrow(out.mat),Nmc,replace=TRUE)
          
          initialXs <- rbeta(prow,dataFinal$x1.a[1],dataFinal$x1.b[1])
          
          print(mean(initialXs))
          
          days <- seq(1,NT)
          ysDet <- forecastStep(IC=mean(initialXs),b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),b4=mean(b4),b5=mean(b5),
                                n=1,NT=length(days),Tair=dataFinal$TairMu[,1],D=dataFinal$D[,1]) #Dependent on IC, but independent for >=1
          ysParam <- forecastStep(IC=mean(initialXs),b0=b0[prow],b1=b1[prow],b2=b2[prow],b3=b3[prow],b4=b4[prow],b5=b5[prow],
                                  n=Nmc,NT=NT,Tair=dataFinal$TairMu[,1],D=dataFinal$D[,1])
          ysProc <- forecastStep(IC=mean(initialXs),b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),b4=mean(b4),b5=mean(b5),
                                 n=Nmc,NT=NT,Q=sqrt(1/out.mat$p.proc[prow]),Tair=dataFinal$TairMu[,1],D=dataFinal$D[,1])
          ysIC <- forecastStep(IC=initialXs,b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),b4=mean(b4),b5=mean(b5),
                               n=length(initialXs),NT=NT,Tair=dataFinal$TairMu[,1],D=dataFinal$D[,1])
          if(dataFinal$N>1){
            for(i in 2:dataFinal$N){
              initialXs <- rbeta(prow,dataFinal$x1.a[i],dataFinal$x1.b[i])
              
              ysDet <- cbind(ysDet,forecastStep(IC=mean(initialXs),b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),b4=mean(b4),b5=mean(b5),
                                                n=1,NT=length(days),Tair=dataFinal$TairMu[,i],D=dataFinal$D[,i]))
              ysParam <- cbind(ysParam,forecastStep(IC=mean(initialXs),b0=b0[prow],b1=b1[prow],b2=b2[prow],b3=b3[prow],b4=b4[prow],b5=b5[prow],
                                                    n=Nmc,NT=NT,Tair=dataFinal$TairMu[,i],D=dataFinal$D[,i]))
              ysProc <- cbind(ysProc,forecastStep(IC=mean(initialXs),b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),b4=mean(b4),b5=mean(b5),
                                                  n=Nmc,NT=NT,Q=sqrt(1/out.mat$p.proc[prow]),Tair=dataFinal$TairMu[,i],D=dataFinal$D[,i]))
              ysIC <- cbind(ysIC,forecastStep(IC=initialXs,b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),b4=mean(b4),b5=mean(b5),
                                              n=length(initialXs),NT=NT,Tair=dataFinal$TairMu[,i],D=dataFinal$D[,i]))
            }
          }
          N.IP.ci = apply(ysParam,2,quantile,c(0.025,0.5,0.975))
          N.Proc.ci = apply(ysProc,2,quantile,c(0.025,0.5,0.975))
          N.IC.ci = apply(ysIC,2,quantile,c(0.025,0.5,0.975))
          predictedValues <- list(ysDet=ysDet,N.IP.ci=N.IP.ci,N.Proc.ci=N.Proc.ci,N.IC.ci=N.IC.ci)
          save(file=outFileName,predictedValues)
        }
      }
    }
  }
}

filePattern <- "_meanTemp_b3_calibration_varBurn.RData"
pdfName <- "UncertaintyAnaylisHindcast_meanTemp_b3.pdf"
uncertaintyAnalysisHindcast(filePattern=filePattern,pdfName=pdfName,b="b3",tID="meanTemp")
filePattern <- "_meanTemp_b4_calibration_varBurn.RData"
pdfName <- "UncertaintyAnaylisHindcast_meanTemp_b4.pdf"
uncertaintyAnalysisHindcast(filePattern=filePattern,pdfName=pdfName,b="b4",tID="meanTemp")
filePattern <- "_dayTemp_b3_calibration_varBurn.RData"
pdfName <- "UncertaintyAnaylisHindcast_dayTemp_b3.pdf"
uncertaintyAnalysisHindcast(filePattern=filePattern,pdfName=pdfName,b="b3",tID="dayTemp")
filePattern <- "_dayTemp_b4_calibration_varBurn.RData"
pdfName <- "UncertaintyAnaylisHindcast_dayTemp_b4.pdf"
uncertaintyAnalysisHindcast(filePattern=filePattern,pdfName=pdfName,b="b4",tID="dayTemp")
