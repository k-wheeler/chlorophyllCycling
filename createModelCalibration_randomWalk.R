##Create Calibration Fits for mean air temperature no trigger model
library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(suncalc)
library(rnoaa)
library(doParallel)

n.cores <- 25
registerDoParallel(cores=n.cores)

dataDirectory <- "data/"
siteData <- read.csv('/projectnb/dietzelab/kiwheel/chlorophyllCycling/allPhenocamDBsitesComplete.csv',header=TRUE)
#siteData <- read.csv('allPhenocamDBsitesComplete.csv',header=TRUE)
sites <- c("harvard","umichbiological","bostoncommon","coweeta","howland2",
           "morganmonroe","missouriozarks","queens","dukehw","lacclair","bbc1","NEON.D08.DELA.DP1.00033",
           "arbutuslake","bartlettir","proctor","oakridge1","hubbardbrook","canadaOA","alligatorriver","readingma",
           "bullshoals","thompsonfarm2N","ashburnham","shalehillsczo")
yearsRemoved <- c(2015,2010,2020,2015,2017,2015,2015,2010,2017,2019,2018,2017,
                  2012,2019,2019,2010,2014,2015,2017,2018,2016,2011,2012,2019)
nchain=5

variableNames <- c("p.PC","x","p.proc")

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
x[i,yr] ~ dnorm(x[(i-1),yr],p.proc)
}
}

#### Priors
for(yr in 1:N){ ##Initial Conditions
x[1,yr] ~ dbeta(x1.a[yr],x1.b[yr]) I(0.001,0.999)
}
p.PC ~ dgamma(s1.PC,s2.PC)
p.proc ~ dgamma(s1.proc,s2.proc)
}
"
ns <- c(49,56,63,70,77,84,91,98,105,112,119,126)#,183)
foreach(s =1:length(sites)) %dopar% {
#for(s in 3:length(sites)){
  siteName <- sites[s]
  print(siteName)
  yearRemoved <- yearsRemoved[s]
  load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData"))
  for(n in ns){
  #outputFileName <- paste0(siteName,"_summer_meanTemp_expBreak_b3_calibration_varBurn.RData")
  outputFileName <- paste0(siteName,"_randomWalk_",n,"_forecast_calibration_varBurn.RData")
  if(!file.exists(outputFileName)){
    #Remove year
    yearInt <- which(dataFinal$years==yearRemoved)
    dataFinal$p[,yearInt] <- NA
    dataFinal$p[(n+1):nrow(dataFinal$p),] <- NA
    # dataFinal$p <- dataFinal$p[,-yearInt]
    # dataFinal$TairMu <- dataFinal$TairMu[,-yearInt]
    # dataFinal$TairPrec <- dataFinal$TairPrec[,-yearInt]
    # dataFinal$TairMuDay <- dataFinal$TairMuDay[,-yearInt]
    # dataFinal$TairPrecDay <- dataFinal$TairPrecDay[,-yearInt]
    # dataFinal$D <- dataFinal$D[,-yearInt]
    # dataFinal$x1.a <- dataFinal$x1.a[-yearInt]
    # dataFinal$x1.b <- dataFinal$x1.b[-yearInt]
    # dataFinal$N <- dataFinal$N - 1
    
    ##Add priors
    dataFinal$s1.PC <- 1.56
    dataFinal$s2.PC <- 0.016
    dataFinal$s1.proc <- 1.56
    dataFinal$s2.proc <- 0.016
    #dataFinal$n <- 46 ##Set for summer
    
    j.model   <- jags.model(file = textConnection(generalModel),
                            data = dataFinal,
                            n.chains = nchain,
                            n.adapt = 1500)
    
    out.burn <- runForecastIter(j.model=j.model,variableNames=variableNames,
                                baseNum = 15000,iterSize = 5000,effSize = 5000,
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
  }
}



