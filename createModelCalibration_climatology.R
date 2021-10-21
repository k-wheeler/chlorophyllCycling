##Create Calibration Fits for mean air temperature no trigger model
library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(suncalc)
library(rnoaa)
library(doParallel)

n.cores <- 12
registerDoParallel(cores=n.cores)

dataDirectory <- "data/"
siteData <- read.csv('/projectnb/dietzelab/kiwheel/chlorophyllCycling/allPhenocamDBsitesComplete.csv',header=TRUE)
#siteData <- read.csv('allPhenocamDBsitesComplete.csv',header=TRUE)
sites <- c("harvard","umichbiological","bostoncommon","coweeta","howland2",
           "morganmonroe","missouriozarks","queens","dukehw","lacclair","bbc1","NEON.D08.DELA.DP1.00033")
yearsRemoved <- c(2015,2010,2020,2015,2017,2015,2015,2010,2017,2019,2018,2017)
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
  siteName <- sites[s]
  print(siteName)
  yearRemoved <- yearsRemoved[s]
  load(paste0(dataDirectory,siteName,"_dataFinal.RData"))
  #outputFileName <- paste0(siteName,"_summer_meanTemp_expBreak_b3_calibration_varBurn.RData")
  outputFileName <- paste0(siteName,"_climatology_calibration_varBurn.RData")
  if(!file.exists(outputFileName)){
    #Remove year
    yearInt <- which(dataFinal$years==yearRemoved)
    dataFinal$p <- dataFinal$p[,-yearInt]
    dataFinal$TairMu <- dataFinal$TairMu[,-yearInt]
    dataFinal$TairPrec <- dataFinal$TairPrec[,-yearInt]
    dataFinal$D <- dataFinal$D[,-yearInt]
    dataFinal$x1.a <- dataFinal$x1.a[-yearInt]
    dataFinal$x1.b <- dataFinal$x1.b[-yearInt]
    dataFinal$N <- dataFinal$N - 1
    
    ##Add priors
    dataFinal$s1.PC <- 1.56
    dataFinal$s2.PC <- 0.016
    #dataFinal$n <- 46 ##Set for summer
    
    j.model   <- jags.model(file = textConnection(generalModel),
                            data = dataFinal,
                            n.chains = nchain,
                            n.adapt = 1500)
    
    out.burn <- runMCMC_Model(j.model=j.model,variableNames=variableNames,
                                baseNum = 15000,iterSize = 5000,effSize = 5000)
    
    ##Thin the data:
    out.mat <- as.matrix(out.burn$params)
    thinAmount <- round(nrow(out.mat)/5000,digits=0)
    out.burn2 <- list()
    out.burn2$params <- window(out.burn$params,thin=thinAmount)
    out.burn2$predict <- window(out.burn$predict,thin=thinAmount)
    out.burn <- out.burn2
    save(out.burn,file = outputFileName)
  }
}



