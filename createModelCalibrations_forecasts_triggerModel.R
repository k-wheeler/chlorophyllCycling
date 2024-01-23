##Create Calibration Fits for mean air temperature no trigger model
library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(suncalc)
library(rnoaa)
library(doParallel)
source('generalVariables.R')

createTriggerModelCalibration <- function(ns,sVals){
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
    #Tair[i,yr] ~ dnorm(TairMu[i,yr],TairPrec[i,yr])
    
    offsetRaw[i,yr] <- max(Tb-TairMu[i,yr],0)
    offset[i,yr] <- ifelse(D[i,yr]<Dstart,offsetRaw[i,yr],0)
    CDD[i,yr] <- CDD[(i-1),yr] + offset[i,yr] * (D[i,yr]/Dstart)
    xmu[i,yr] <- max(min((ifelse(CDD[i,yr]>CDDcrit,x[(i-1),yr] - summerRate, x[(i-1),yr] - fallRate)),x[1,yr]),0)

    x[i,yr] ~ dnorm(xmu[i,yr],p.proc)
    }
    }
    
    #### Priors
    for(yr in 1:N){ ##Initial Conditions
    x[1,yr] ~ dbeta(x1.a[yr],x1.b[yr]) I(0.001,0.999)
    CDD[1,yr] <- 0
    }
    p.PC ~ dgamma(s1.PC,s2.PC)
    p.proc ~ dgamma(s1.proc,s2.proc)
    
    Dstart <- 15
    CDDcrit ~ dunif(0,2000)
    Tb <- 20
    summerRate <- 0.002 #~ dunif(0,1)
    fallRate ~ dunif(0,1)
    
  }
    "

  registerDoParallel(cores=n.cores)
  
  variableNames <- c("p.PC","x","p.proc","CDDcrit","fallRate")
  
  foreach(s =sVals) %dopar% {
    siteName <- sites[s]
    print(siteName)
    yearRemoved <- yearsRemoved[s]
    load(paste0(dataDirectory,siteName,"_dataFinal.RData"))
    
    for(n in ns){
      print(n)
      outputFileName <- paste0(siteName,"_",n,"basicTrigger_calibration_varBurn.RData")
      
      dataFinal$p[(n+1):nrow(dataFinal$p),] <- NA
      
      yearInt <- which(dataFinal$years==yearRemoved)
      dataFinal$p[,yearInt] <- NA
      
      ##Add priors
      dataFinal$s1.PC <- 1.56
      dataFinal$s2.PC <- 0.016
      dataFinal$s1.proc <- 1.56
      dataFinal$s2.proc <- 0.016
      
      #save(file=initsFileName,inits) #Need to save inits for dic calculations
      j.model <- try(jags.model(file = textConnection(generalModel),
                                data = dataFinal,
                                n.chains = nchain,
                                n.adapt = 1000))#Load Model Output 
      if(inherits(j.model,"try-error")){
        next
      }
      
      out.burn <- try(runForecastIter(j.model=j.model,variableNames=variableNames,
                                      baseNum = 3000,iterSize = 1000,effSize = 1000, maxIter=1000000,
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


createTriggerModelCalibration(ns=182,sVals=1) #Change for number of included days (ns) and sites (sVals)
