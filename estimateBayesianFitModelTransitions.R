library('rjags')
library('runjags')
library('ecoforecastR')
library('PhenologyBayesModeling')
library(doParallel)

n.cores <- 4
registerDoParallel(cores=n.cores)

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
#' @param D
#' @import rjags
#' @import runjags
#' @import ecoforecastR
#'
#' @return
#' @export
#'
#' @examples
forecastStep <- function(IC,b0,b1,b2,b3,Q=0,n,NT,Tair,D){
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

##' Create a Bayes Model for a deciduous broadleaf site
##'
##' @param yobs
##' @import rjags
##' @import runjags
##' @import PhenologyBayesModeling
##' @export
createBayesModel_Fall <- function(yobs) {
  nchain = 5
  data <- list(yobs=yobs,x=seq(1,153),n=153)
  
  data$c.alp <- 8
  data$c.bet <- 1.5
  data$d <- 0
  data$s1 <- 0.001
  data$s2 <- 0.00001
  data$p.Tran <- 1/(40**2)
  data$p.b <- 1/(0.05**2)
  data$mean.TranF <- 300
  data$mean.bF <- 0.10
  data$obs.prec <- 100
  
  DB_model <- "
    model{
    ##priors
    TranF ~ dnorm(mean.TranF,p.Tran)  ##F for fall/autumn
    bF ~ dnorm(mean.bF,p.b)
    prec ~ dgamma(s1,s2)
    c ~ dbeta(c.alp,c.bet)
    for(i in 1:n){
    mu[i] <- c/(1+exp(bF*(x[i]-TranF)))+d ##process model for fall
    y[i] ~ dnorm(mu[i],prec)
    yobs[i] ~ dnorm(y[i],obs.prec)
    }
    }
    "
  
  j.model   <- jags.model(file = textConnection(DB_model),
                          data = data,
                          n.chains = nchain)
  return(j.model)
}

load('harvard_b3_final_calibration_varBurn.RData')

Nmc <- 1000
out.mat <- data.frame(as.matrix(out.burn$param))

b0 <- out.mat$b0
b1 <- out.mat$b1
b2 <- out.mat$b2
b3 <- out.mat$b3
prow = sample.int(nrow(out.mat),Nmc,replace=TRUE)

#files <- dir(pattern = "_b3_final_calibration_varBurn.RData")

sites <- c("asuhighlands","bullshoals","macleish","harvard")

minCt <- 1

#for(f in 1:length(sites)){
foreach(f=1:length(sites)) %dopar% {
  fileName <- paste0(sites[f],"_b3_final_calibration_varBurn.RData")
  print(fileName)
  tle <- strsplit(fileName,"_")[[1]][1:2]
  siteName <- strsplit(fileName,"_")[[1]][1]
  
  baseTempOrig <- 20
  load(fileName)
  load(paste(siteName,"_CDD",baseTempOrig,"_dataFinal.RData",sep=""))
  
  if(ncol(dataFinal$TairMu)>1){
    out.mat <- data.frame(as.matrix(out.burn$param))
    NT <- length(dataFinal$TairMu[,1])
    out.mat.pred <- data.frame(as.matrix(out.burn$predict))
    
    prow.pred = sample.int(nrow(out.mat),Nmc,replace=TRUE)
    initialXs <- out.mat.pred[prow.pred,1]
    days <- seq(1,NT)
    
    for(i in 1:min(dataFinal$N,4)){
      initialXs <- out.mat.pred[prow.pred,dataFinal$n*(i-1)+1]
      ysDet <- as.numeric(forecastStep(IC=mean(initialXs),b0=mean(b0),b1=mean(b1),b2=mean(b2),
                                       b3=mean(b3),n=1,NT=length(days),Tair=dataFinal$TairMu[,i],
                                       D=dataFinal$D[,1]))
      j.model <- createBayesModel_Fall(yobs=ysDet)
      var.burn <- runMCMC_Model(j.model = j.model,variableNames = c("TranF","bF","prec","c"),baseNum=100000,
                                iterSize = 10000,sampleCutoff = 3000)
      save(var.burn,file=paste0(siteName,"_",i,"_ysDet_phenoCurve_varBurn.RData"))
    }
  }
}
