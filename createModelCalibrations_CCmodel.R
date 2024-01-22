##Create Calibration Fits for mean air temperature no trigger model
library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(suncalc)
library(rnoaa)
library(doParallel)
source('generalVariables.R')

createChlorophyllCyclingModelCalibration <- function(summerOnly=FALSE,ns,sVals){
  print(ns)
  generalModel = "
model {
    ### Data Models for complete years
    for(yr in 1:(N)){
    for(i in 1:n){
    p[i,yr] ~ dnorm(x[i,yr],p.PC)
    }
    }
    
    #### Process Model
    for(yr in 1:(N)){
    for(i in 2:n){
    Tair[i,yr] ~ dnorm(TairMu[i,yr],TairPrec[i,yr])
    
    xmu[i,yr] <- max(min((x[(i-1),yr] + b4 * x[(i-1),yr]) + max(0,(b0 + b3 * Tair[i,yr] * D[i,yr])),x[1,yr]),0)
    x[i,yr] ~ dnorm(xmu[i,yr],p.proc)
    }
    }
    
    #### Priors
    for(yr in 1:N){ ##Initial Conditions
    x[1,yr] ~ dbeta(x1.a[yr],x1.b[yr]) I(0.001,0.999)
    }
    p.PC ~ dgamma(s1.PC,s2.PC)
    p.proc ~ dgamma(s1.proc,s2.proc)
    
    b0 ~ dunif(b0_lower,b0_upper)
    b3 ~ dunif(b3_lower,b3_upper)
    b4 ~ dunif(b4_lower,b4_upper)
    
  }"
  
  registerDoParallel(cores=n.cores)
  variableNames <- c("p.PC","x","p.proc","b0","b3","b4")
  
  foreach(s =sVals) %dopar% {
    siteName <- sites[s]
    print(siteName)
    yearRemoved <- yearsRemoved[s]
    load(paste0(dataDirectory,siteName,"_dataFinal.RData"))
    
    for(n in ns){
      print(n)
      outputFileName <- paste0(CCmodelOutputsFolder,siteName,"_",n,"_ccModel_forecast_calibration_varBurn.RData")
      dataFinal$p[(n+1):nrow(dataFinal$p),] <- NA

      #Remove year
      yearInt <- which(dataFinal$years==yearRemoved)
      dataFinal$p[,yearInt] <- NA
      
      ##Add priors
      dataFinal$s1.PC <- 1.56
      dataFinal$s2.PC <- 0.016
      dataFinal$s1.proc <- 1.56
      dataFinal$s2.proc <- 0.016
      dataFinal$b0_lower <- -1 #intercept of synthesis
      dataFinal$b0_upper <- 1 #intercept of synthesis
      
      dataFinal$b3_lower <- 0 #slope of synthesis
      dataFinal$b3_upper <- 1
      dataFinal$b4_lower <- -1 #Breakdown rate
      dataFinal$b4_upper <- 0 #Breakdown rate
      
      j.model <- try(jags.model(file = textConnection(generalModel),
                                data = dataFinal,
                                n.chains = nchain,
                                n.adapt = 3000))#Load Model Output 
      if(inherits(j.model,"try-error")){
        next
      }
      
      out.burn <- try(runForecastIter(j.model=j.model,variableNames=variableNames,
                                      baseNum = 15000,iterSize = 5000,effSize = 5000, maxIter=1000000,
                                      partialFile = paste("partial_",outputFileName,sep="")))
      if(inherits(out.burn,"try-error")){
        next
      }
      
      ##Thin the data:
      if(typeof(out.burn)!=typeof(FALSE)){
        out.mat <- as.matrix(out.burn$params)
        thinAmount <- round(nrow(out.mat)/5000,digits=0)
        out.burn2 <- list()
        out.burn2$params <- window(out.burn$params,thin=thinAmount)
        out.burn2$predict <- window(out.burn$predict,thin=thinAmount)
        out.burn <- out.burn2
      }
      save(out.burn,file = outputFileName)
      print(paste("saved:",outputFileName))
    }
  }
  
}

createChlorophyllCyclingModelCalibration(ns=seq(45,183),
                                         sVals=seq_along(sites)) #Change for number of included days (ns) and sites (sVals)
