##Create Calibration Fits for mean air temperature no trigger model
library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(suncalc)
library(rnoaa)
library(doParallel)
source('generalVariables.R')

registerDoParallel(cores=n.cores)

tranOffsets <- read.csv(allPhenoTranFile,header=TRUE)

variableNames <- c("p.PC","mu")

generalModel = "
model {
### Data Models for complete years

for(i in 1:n){
  mu[i] ~ dunif(0,1) #prior on means 
  for(yr in 1:N){
    p[i,yr] ~ dnorm(mu[i],p.PC) #data model 
  }
}
p.PC ~ dgamma(s1.PC,s2.PC)

}
"
n=184

foreach(s =1:length(sites)) %dopar% {
  siteName <- as.character(sites[s])
  print(siteName)
  
  load(paste0(dataDirectory,siteName,"_dataFinal.RData"))
  outputFileName <- paste0(climatologyModelOutputsFolder,siteName,"_climatology_forecast_calibration_varBurn.RData")
  
  if(!file.exists(outputFileName)){
    #Remove year
    if(!is.na(dataFinal$yearRemoved)){
      yearRemoved <- dataFinal$yearRemoved
      yearInt <- which(dataFinal$years==yearRemoved)
      
      dataFinal$p[,yearInt] <- NA
      dataFinal$n <- n
    }
    
    ##Add priors
    dataFinal$s1.PC <- 1.56
    dataFinal$s2.PC <- 0.016
    
    j.model   <- jags.model(file = textConnection(generalModel),
                            data = dataFinal,
                            n.chains = nchain,
                            n.adapt = 1500)
    
    out.burn <- runMCMC_Model(j.model=j.model,variableNames=variableNames,
                              baseNum = 5000,iterSize = 2000,maxGBR=50)
    
    ##Thin the data:
    out.mat <- as.matrix(out.burn)
    thinAmount <- round(nrow(out.mat)/5000,digits=0)
    out.burn2 <- list()
    out.burn2 <- window(out.burn,thin=thinAmount)
    out.burn <- out.burn2
    save(out.burn,file = outputFileName)
  }
}



