library('rjags')
library('runjags')
library('ecoforecastR')
library('PhenologyBayesModeling')
library(doParallel)

n.cores <- 25
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
forecastStep <- function(IC,b0,b1,b2,b3,b4,Q=0,n,NT,Tair,D){
  x <- matrix(NA,n,NT)
  if(length(IC)==1){
    Xprev <- rep(IC,n)
  }else{
    Xprev <- IC
  }
  
  for(t in 1:NT){
    bd <- Xprev + b4 * Xprev #b0 + (b1 * Xprev) + (b2 * Xprev **2) 
    syn <- b0 + b1 * Tair[t] + b2 * D[t] + b3 * Tair[t] * D[t] 
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
createChangepointModel_Fall <- function(yobs) {
  nchain = 5
  data <- list(yobs=yobs,x=seq(1,length(yobs)),n=length(yobs))
  data$obs.prec <- 100
  data$s1 <- 0.001
  data$s2 <- 0.00001
  data$c.alp <- 8
  data$c.bet <- 1.5
  DB_model <- "
  model{
  ##priors
  #prec ~ dgamma(s1,s2)
  y[1] ~ dbeta(c.alp,c.bet)
  mS ~ dunif(0,1)
  mF ~ dunif(0,1)
  k ~ dunif(0,100)

  for(i in 2:n){
  muS[i] <- y[(i-1)] - mS
  muF[i] <- y[(i-1)] - mF
  y[i] <- ifelse(x[i]>k,muF[i],muS[i])
  #y[i] ~ dnorm(mu[i],prec)
  yobs[i] ~ dnorm(y[i],obs.prec)
  }
  }
  "
  
  j.model   <- jags.model(file = textConnection(DB_model),
                          data = data,
                          n.chains = nchain)
  return(j.model)
}
runSites <- c("harvard","umichbiological","bostoncommon","coweeta","howland2",
           "morganmonroe","missouriozarks","queens","dukehw","lacclair","bbc1","NEON.D08.DELA.DP1.00033",
           "arbutuslake","bartlettir","proctor","oakridge1","hubbardbrook","asa","canadaOA","alligatorriver","readingma",
           "bullshoals","thompsonfarm2N","ashburnham","shalehillsczo")
site <- 1
load(paste0(runSites[site],'_meanTemp_fall_expBreak_slope_b3_calibration_varBurn.RData'))

Nmc <- 1000
out.mat <- data.frame(as.matrix(out.burn$param))
b <- "b3"
b0 <- out.mat$b0
b4 <- out.mat$b4
if(b=="b1"){
  b1 <- out.mat$b1
  b2 <- rep(0,length(b1))
  b3 <- rep(0,length(b1))
}else if(b=="b2"){
  b2 <- out.mat$b2
  b1 <- rep(0,length(b2))
  b3 <- rep(0,length(b2))
}else if(b=="b3"){
  b3 <- out.mat$b3
  b1 <- rep(0,length(b3))
  b2 <- rep(0,length(b3))
}
prow = sample.int(nrow(out.mat),Nmc,replace=TRUE)

siteData <- read.csv('/projectnb/dietzelab/kiwheel/chlorophyllCycling/allPhenocamDBsitesComplete.csv',header=TRUE)

sites <- as.character(siteData$siteName)
dataDirectory <- "data/"

#for(f in 1:length(sites)){
foreach(s=1:length(sites)) %dopar% {
  siteName <- sites[s]
  print(siteName)
  
  load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData")) #Load Data
  
  if(ncol(dataFinal$TairMu)>1){
    NT <- length(dataFinal$TairMu[,1])

    #initialXs <- rbeta(prow,dataFinal$x1.a[1],dataFinal$x1.b[1])
    days <- seq(1,NT)
    
    for(i in 1:dataFinal$N){
      if(!file.exists(paste0(siteName,"_",dataFinal$years[i],"_ysDet_changePointCurve_",runSites[site],"_meanTemp_fall_b1_varBurn.RData")))
      initialXs <- rbeta(prow,dataFinal$x1.a[i],dataFinal$x1.b[i])
      ysDet <- as.numeric(forecastStep(IC=mean(initialXs),b0=mean(b0),b1=mean(b1),b2=mean(b2),
                                       b3=mean(b3),b4=mean(b4),n=1,NT=length(days),Tair=dataFinal$TairMu[,i],
                                       D=dataFinal$D[,i]))
      if(length(which(ysDet<0.10))>0){
        ysDet <- ysDet[1:(which(ysDet<0.10)[1]-1)]
      }
      j.model <- createChangepointModel_Fall(yobs=ysDet)
      variables <- c("mS","mF","y[1]","k")
      var.burn <- runMCMC_Model(j.model = j.model,variableNames = variables, baseNum=20000,
                                iterSize = 10000,sampleCutoff = 2000)
      save(var.burn,file=paste0(siteName,"_",dataFinal$years[i],"_ysDet_changePointCurve_",runSites[site],"_meanTemp_fall_b1_varBurn.RData"))
    }
  }
}
