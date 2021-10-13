##Create Calibration Fits for mean day air temperature no trigger model
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

variableNames <- c("p.PC","x","b0","b1","b2","b5","p.proc")

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

xmu[i,yr] <- max(min(x[(i-1),yr] + (b0 + (b1 * x[(i-1),yr]) + (b2 * x[(i-1),yr] ** 2)) + max(0,(b5 * D[i,yr])),x[1,yr]),0)
x[i,yr] ~ dnorm(xmu[i,yr],p.proc)
}
}

#### Priors
for(yr in 1:N){ ##Initial Conditions
x[1,yr] ~ dbeta(x1.a[yr],x1.b[yr])
}
p.PC ~ dgamma(s1.PC,s2.PC)
p.proc ~ dgamma(s1.proc,s2.proc)

b0 ~ dunif(b0_lower,b0_upper)
b1 ~ dunif(b1_lower,b1_upper)
b2 ~ dunif(b2_lower,b2_upper)
b5 ~ dunif(b5_lower,b5_upper)
}
"

foreach(s =1:length(sites)) %dopar% {
  siteName <- sites[s]
  print(siteName)
  yearRemoved <- yearsRemoved[s]
  load(paste0(dataDirectory,siteName,"_dataFinal.RData"))
  outputFileName <- paste0(siteName,"_summer_b5_calibration_varBurn.RData")
  
  #Remove year
  yearInt <- which(dataFinal$years==yearRemoved)
  dataFinal$p <- dataFinal$p[,-yearInt]
  dataFinal$TairMuDay <- dataFinal$TairMuDay[,-yearInt]
  dataFinal$TairPrecDay <- dataFinal$TairPrecDay[,-yearInt]
  dataFinal$D <- dataFinal$D[,-yearInt]
  dataFinal$x1.a <- dataFinal$x1.a[-yearInt]
  dataFinal$x1.b <- dataFinal$x1.b[-yearInt]
  dataFinal$N <- dataFinal$N - 1
  
  ##Add priors
  dataFinal$s1.PC <- 1.56
  dataFinal$s2.PC <- 0.016
  dataFinal$s1.proc <- 1.56
  dataFinal$s2.proc <- 0.016
  dataFinal$b0_lower <- -1
  dataFinal$b0_upper <- 0
  dataFinal$b1_lower <- -1
  dataFinal$b1_upper <- 0
  dataFinal$b2_lower <- -1
  dataFinal$b2_upper <- 0
  dataFinal$b5_lower <- 0
  dataFinal$b5_upper <- 1
  dataFinal$n <- 46 #set for summer 
  
  j.model   <- jags.model(file = textConnection(generalModel),
                          data = dataFinal,
                          n.chains = nchain,
                          n.adapt = 2000)
  
  out.burn <- runForecastIter(j.model=j.model,variableNames=variableNames,
                              baseNum = 40000,iterSize = 10000,effSize = 5000,
                              partialFile = paste("partial_",outputFileName,sep=""))

  ##Thin the data:
  out.mat <- as.matrix(out.burn$params)
  thinAmount <- round(nrow(out.mat)/5000,digits=0)
  out.burn2 <- list()
  out.burn2$params <- window(out.burn$params,thin=thinAmount)
  out.burn2$predict <- window(out.burn$predict,thin=thinAmount)
  out.burn <- out.burn2
  save(out.burn,file = outputFileName)
}



