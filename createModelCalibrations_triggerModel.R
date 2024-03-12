##Create Calibration Fits for mean air temperature no trigger model
library(rjags)
library(runjags)
library(doParallel)
source('generalVariables.R')
source('runModelIterations.R')

allTrans <- read.csv(allPhenoTranFile)

createTriggerModelCalibration <- function(sVals){
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
    
    offsetRaw[i,yr] <- max(Tb-Tair[i,yr],0)
    offset[i,yr] <- ifelse(D[i,yr]<Dstart,offsetRaw[i,yr],0)
    CDD[i,yr] <- CDD[(i-1),yr] + offset[i,yr] * (D[i,yr]/Dstart)
    #xmu[i,yr] <- max(min((ifelse(CDD[i,yr]>CDDcrit,x[(i-1),yr] - summerRate, x[(i-1),yr] - fallRate)),x[1,yr]),0)

    xmu[i,yr] <- max(min((ifelse(CDD[i,yr]>CDDcrit,x[(i-1),yr] - summerRate, x[(i-1),yr] + fallRate * x[(i-1),yr] * (1-x[(i-1),yr]))),x[1,yr]),0)   
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
    
    Dstart ~ dunif(11,17)
    CDDcrit ~ dunif(0,500)
    Tb ~ dunif(5,25) #<- 20
    summerRate ~ dunif(0,0.1)
    fallRate ~ dunif(-0.99,0)
    
  }
    "
  registerDoParallel(cores=min(n.cores,length(sVals)))
  
  variableNames <- c("p.PC","x","p.proc","CDDcrit","fallRate","summerRate","Dstart","Tb")
  
  foreach(s =sVals) %dopar% {
    siteName <- sites[s]
    print(siteName)
    yearRemoved <- yearsRemoved[s]
    load(paste0(dataDirectory,siteName,"_dataFinal.RData"))
    
    outputFileName <- paste0(triggerModelOutputsFolder,siteName,"_triggerModel_calibration_varBurn.RData")
    partialFileName <- paste0(triggerModelOutputsFolder,siteName,"_triggerModel_calibration_varBurn_partial.RData")
    
    if(!file.exists(outputFileName)){
      approximateTran <- round(allTrans[allTrans$siteName==siteName,'meanDOY']-181,digits=0)
      
      mdl <- summary(lm(dataFinal$p[approximateTran:183,1]~seq(approximateTran,183)))
      fallRateInit <- abs(mdl$coefficients[2,1])
      
      mdl <- summary(lm(dataFinal$p[1:approximateTran,1]~seq(1,approximateTran)))
      summerRateInit <- abs(mdl$coefficients[2,1])
      
      yearInt <- which(dataFinal$years==yearRemoved)
      dataFinal$p[,yearInt] <- NA
      
      ##Add priors
      dataFinal$s1.PC <- 1.56
      dataFinal$s2.PC <- 0.016
      dataFinal$s1.proc <- 1.56
      dataFinal$s2.proc <- 0.016
      
      inits <- list()
      
      if(siteName=="harvard"){
        fallRateInit <- -0.2
        TbInit <- 15
        CDDcritInit <- 25
        DstartInit <- 16
      }else if(siteName=="coweeta"){
        fallRateInit <- -0.2
        TbInit <- 18
        CDDcritInit <- 50
        DstartInit <- 15
      }else if(siteName=="missouriozarks"){
        fallRateInit <- -0.15
        TbInit <- 18
        CDDcritInit <- 100
        DstartInit <- 15
      }else if(siteName=="NEON.D08.DELA.DP1.00033"){
        fallRateInit <- -0.05
        TbInit <- 23
        CDDcritInit <- 15
        DstartInit <- 16
      }else if(siteName=="oakridge1"){
        fallRateInit <- -0.15
        TbInit <- 20
        CDDcritInit <- 50
        DstartInit <- 16
      }else if(siteName=="dukehw"){
        fallRateInit <- -0.08
        TbInit <- 20
        CDDcritInit <- 25
        DstartInit <- 16
      }else if(siteName=="umichbiological"){
        fallRateInit <- -0.1
        TbInit <- 15
        CDDcritInit <- 50
        DstartInit <- 16
        summerRateInit <- 0.001
      }else if(siteName=="proctor"){
        fallRateInit <- -0.17
        TbInit <- 15
        CDDcritInit <- 50
        DstartInit <- 16
      }else if(siteName=="hubbardbrook"){
        fallRateInit <- -0.17
        TbInit <- 18
        CDDcritInit <- 50
        DstartInit <- 16
      }else if(siteName=="howland2"){
        fallRateInit <- -0.17
        TbInit <- 18
        CDDcritInit <- 50
        DstartInit <- 16
      }else{
        fallRateInit <- -0.2
        TbInit <- 15
        CDDcritInit <- 50
        DstartInit <- 16
      }
      for(c in 1:nchain){
        inits[[c]] <- list(summerRate=rnorm(1,summerRateInit,0.0001),
                           fallRate=rnorm(1,fallRateInit,0.0001),
                           Tb=min(rnorm(1,TbInit,0.5),19.9),
                           CDDcrit=rnorm(1,CDDcritInit,5),
                           Dstart=min(rnorm(1,DstartInit,0.5),16.5))
      }
      
      #save(file=initsFileName,inits) #Need to save inits for dic calculations
      j.model <- try(jags.model(file = textConnection(generalModel),
                                data = dataFinal,
                                n.chains = nchain,
                                inits = inits,
                                n.adapt = 2000))#Load Model Output 
      if(inherits(j.model,"try-error")){
        next
      }
      
      out.burn <- try(runForecastIter(j.model=j.model,variableNames=variableNames,
                                      baseNum = 10000,iterSize = 5000,effSize = 5000, maxIter=1000000,partialFile = partialFileName))
      #partialFile = paste("partial_",outputFileName,sep="")))
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
        save(out.burn,file = outputFileName)
        print(paste("saved:",outputFileName))
      }else{
        print(paste(siteName,"Did not converge"))
      }

    }
  }
}

#c(1,11,17,12,20)
createTriggerModelCalibration(sVals=seq_along(sites)) #Change for number of included days (ns) and sites (sVals)
#c(9,7,12,16,2)
#Playing around
# plot(dataFinal$p[,1],pch=20)
# plot(dataFinal$D[,1],pch=20)
# abline(h=15,col="red")
# 
# Tb <- 20
# offsetRaw <- matrix(nrow=dataFinal$n,ncol=dataFinal$N)
# offset <- matrix(nrow=dataFinal$n,ncol=dataFinal$N)
# CDD <- matrix(nrow=dataFinal$n,ncol=dataFinal$N)
# offsetRaw <- matrix(nrow=dataFinal$n,ncol=dataFinal$N)
# CDD[1,] <- 0
# Tair <- dataFinal$TairMu
# D <- dataFinal$D
# Dstart <- 15
# 
# for(yr in 1:dataFinal$N){
#   for(i in 2:dataFinal$n){
#     offsetRaw[i,yr] <- max(Tb-Tair[i,yr],0)
#     if(D[i,yr]<Dstart){
#       offset[i,yr] <- offsetRaw[i,yr]
#     }else{
#       offset[i,yr] <- 0
#     }
#     
#     CDD[i,yr] <- CDD[(i-1),yr] + offset[i,yr] * (D[i,yr]/Dstart)
#   }
# }
# 
# par(mfrow=c(1,1))
# yr=9
# plot(dataFinal$p[,yr],pch=20)
# abline(mdl,col="red")
# plot(CDD[,yr],pch=20,ylim=c(0,200))
# 
# mdl <- lm(dataFinal$p[1:100,1]~seq(1,100))
# 
# out <- coda.samples(j.model,variable.names = variableNames,n.iter=3000)
