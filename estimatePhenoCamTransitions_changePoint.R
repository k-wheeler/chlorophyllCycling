library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(doParallel)

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

n.cores <- 28

#register the cores.
registerDoParallel(cores=n.cores)

dataDirectory <- "data/"
siteData <- read.csv('/projectnb/dietzelab/kiwheel/chlorophyllCycling/allPhenocamDBsitesComplete.csv',header=TRUE)
siteName <- "harvardblo"
#for(s in 35:nrow(siteData)){
foreach(s=1:nrow(siteData)) %dopar% {
  siteName <- as.character(siteData$siteName[s])
  print(siteName)
  
  load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData"))
  
  for(yr in 1:dataFinal$N){
    yrName <- dataFinal$years[yr]
    if(!file.exists(paste0(siteName,"_",yrName,"_PhenoCam_changePointCurve_varBurn.RData"))){
      #p <- dataFinal$p[dataFinal$p[,yr]>0.15,yr]
      if(length(which(dataFinal$p[,yr]<0.10))>0){
        p <- dataFinal$p[1:(which(dataFinal$p[,yr]<0.10)[1]-1),yr]
      }else{
        p <- dataFinal$p[,yr]
      }

      j.model <- createChangepointModel_Fall(yobs=p)
      variables <- c("mS","mF","y[1]","k")
      var.burn <- runMCMC_Model(j.model = j.model,variableNames = variables, baseNum=100000,
                                iterSize = 50000,sampleCutoff = 5000,maxGBR = 50)
      save(var.burn,file=paste0(siteName,"_",yrName,"_PhenoCam_changePointCurve_varBurn.RData"))
    }else{
      load(paste0(siteName,"_",yrName,"_PhenoCam_changePointCurve_varBurn.RData"))
      if(typeof(var.burn)==typeof(FALSE)){
        if(length(which(dataFinal$p[,yr]<0.10))>0){
          p <- dataFinal$p[1:(which(dataFinal$p[,yr]<0.10)[1]-1),yr]
        }else{
          p <- dataFinal$p[,yr]
        }
        
        j.model <- createChangepointModel_Fall(yobs=p)
        variables <- c("mS","mF","y[1]","k")
        var.burn <- runMCMC_Model(j.model = j.model,variableNames = variables, baseNum=100000,
                                  iterSize = 50000,sampleCutoff = 5000,maxGBR = 50)
        save(var.burn,file=paste0(siteName,"_",yrName,"_PhenoCam_changePointCurve_varBurn.RData"))
      }
    }
  }
}
