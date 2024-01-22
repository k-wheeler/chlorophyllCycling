library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(doParallel)
source('generalVariables.R')

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
  k ~ dunif(0,183)
  
  for(i in 2:n){
  muS[i] <- y[(i-1)] - mS
  muF[i] <- y[(i-1)] - mF
  y[i] <- ifelse(x[i]>k,muF[i],muS[i])
  #y[i] ~ dnorm(mu[i],prec)
  yobs[i] ~ dnorm(y[i],obs.prec)
  }
  }
  "
  inits <- list()
  c("mS","mF","y[1]","k")
  for(i in 1:nchain){
    inits[[i]] <- list(mS = rnorm(1,0.0025,0.0001),
                       mF = rnorm(1,0.03,0.001),
                       k = rnorm(1,140,1))
  }
  
  j.model   <- jags.model(file = textConnection(DB_model),
                          data = data,
                          inits = inits,
                          n.chains = nchain,
                          n.adapt = 2000)
  return(j.model)
}

#register the cores.
registerDoParallel(cores=n.cores)

foreach(s=1:nrow(siteData):1) %dopar% {
  siteName <- as.character(siteData$siteName[s])
  print(siteName)
  
  load(paste0(dataDirectory,siteName,"_dataFinal.RData"))
  
  for(yr in 1:dataFinal$N){
    yrName <- dataFinal$years[yr]
    print(yrName)
    if(!file.exists(paste0(transitionEstimateOutputsFolder,siteName,"_",yrName,"_PhenoCam_changePointCurve_varBurn.RData"))){
      #p <- dataFinal$p[dataFinal$p[,yr]>0.15,yr]
      if(length(which(dataFinal$p[,yr]<0.15))>0){
        p <- dataFinal$p[1:(which(dataFinal$p[,yr]<0.15)[1]-1),yr]
      }else{
        p <- dataFinal$p[,yr]
      }

      j.model <- createChangepointModel_Fall(yobs=p)
      variables <- c("mS","mF","y[1]","k")
      var.burn <- runMCMC_Model(j.model = j.model,variableNames = variables, baseNum=50000,
                                iterSize = 10000,sampleCutoff = 5000,maxGBR = 50)
      if(typeof(var.burn)!=typeof(FALSE)){
        out.mat <- as.matrix(var.burn)
        thinAmount <- round(nrow(out.mat)/5000,digits=0)
        var.burn <- window(var.burn,thin=thinAmount)
      }
      save(var.burn,file=paste0(transitionEstimateOutputsFolder,siteName,"_",yrName,"_PhenoCam_changePointCurve_varBurn.RData"))
    }else{
      load((paste0(transitionEstimateOutputsFolder,siteName,"_",yrName,"_PhenoCam_changePointCurve_varBurn.RData")))
      if(typeof(var.burn)==typeof(FALSE)){
        if(length(which(dataFinal$p[,yr]<0.10))>0){
          p <- dataFinal$p[1:(which(dataFinal$p[,yr]<0.10)[1]-1),yr]
        }else{
          p <- dataFinal$p[,yr]
        }
        
        j.model <- createChangepointModel_Fall(yobs=p)
        variables <- c("mS","mF","y[1]","k")
        var.burn <- runMCMC_Model(j.model = j.model,variableNames = variables, baseNum=100000,
                                  iterSize = 5000,sampleCutoff = 1000,maxGBR = 5)
        if(typeof(var.burn)!=typeof(FALSE)){
          out.mat <- as.matrix(var.burn)
          thinAmount <- round(nrow(out.mat)/5000,digits=0)
          var.burn <- window(var.burn,thin=thinAmount)
        }
        save(var.burn,file=(paste0(transitionEstimateOutputsFolder,siteName,"_",yrName,"_PhenoCam_changePointCurve_varBurn.RData")))
      }
    }
  }
}
