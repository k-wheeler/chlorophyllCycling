library('rjags')
library('runjags')
library('ecoforecastR')
library('PhenologyBayesModeling')
library(doParallel)

n.cores <- 4
registerDoParallel(cores=n.cores)

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

sites <- c("asuhighlands","bullshoals","macleish","harvard")

#for(f in 1:length(sites)){
foreach(f=1:length(sites)) %dopar% {
  fileName <- paste0(sites[f],"_b3_final_calibration_varBurn.RData")
  print(fileName)
  tle <- strsplit(fileName,"_")[[1]][1:2]
  siteName <- strsplit(fileName,"_")[[1]][1]
  
  baseTempOrig <- 20

  load(paste(siteName,"_CDD",baseTempOrig,"_dataFinal.RData",sep=""))
  
  if(ncol(dataFinal$TairMu)>1){

    for(i in 1:min(dataFinal$N,4)){

      j.model <- createBayesModel_Fall(yobs=dataFinal$p[,i])
      variables <- c("TranF","bF","prec","c")
      var.burn <- runMCMC_Model(j.model = j.model,variableNames = variables, baseNum=10000,
                                iterSize = 10000,sampleCutoff = 3000)
      save(var.burn,file=paste0(siteName,"_",i,"_p_phenoCurve_varBurn.RData"))
      #save(var.burn,file=paste0(siteName,"_",i,"_ysDet_changePointCurve_varBurn.RData"))
    }
  }
}