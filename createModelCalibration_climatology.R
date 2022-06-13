##Create Calibration Fits for mean air temperature no trigger model
library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(suncalc)
library(rnoaa)
library(doParallel)

n.cores <- 16
registerDoParallel(cores=n.cores)
summerOnly <- TRUE

offset <- 0 #24
tranOffsets <- read.csv('phenocamTransitions_fromMean.csv',header=TRUE)

dataDirectory <- "data/"
siteData <- read.csv('/projectnb/dietzelab/kiwheel/chlorophyllCycling/allPhenocamDBsitesComplete.csv',header=TRUE)
#siteData <- read.csv('allPhenocamDBsitesComplete.csv',header=TRUE)
sites <- c("harvard","umichbiological","bostoncommon","coweeta","howland2",
           "morganmonroe","missouriozarks","queens","dukehw","lacclair","bbc1","NEON.D08.DELA.DP1.00033",
           "bartlettir","proctor","oakridge1","hubbardbrook","alligatorriver","readingma","bullshoals",
           "willowcreek","downerwoods","laurentides","russellsage","sanford","boundarywaters")
nchain=5

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

foreach(s =1:length(sites)) %dopar% {
#for(s in nrow(siteData):1){
  #siteName <- as.character(siteData$siteName[s])
  siteName <- as.character(sites[s])
  print(siteName)
  
  load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData"))
  if(summerOnly){
  tranID <- which(tranOffsets$siteName==siteName)
  tran_DOY <- as.numeric(tranOffsets[tranID,2])
  
  n=round(tran_DOY+offset,digits=0)-182
  outputFileName <- paste0(siteName,"_climatology_forecast_calibration_beforeTransition_varBurn2.RData")
  }else{
    outputFileName <- paste0(siteName,"_climatology_forecast_calibration_varBurn2.RData")
    n=184
  }
  
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
    #dataFinal$n <- 46 ##Set for summer
    
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



