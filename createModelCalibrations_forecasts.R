##Create Calibration Fits for mean air temperature no trigger model
library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(suncalc)
library(rnoaa)
library(doParallel)

createChlorophyllCyclingModelCalibration <- function(b0=c(-1,0),b1=c(0,0),b2=c(0,0),b3=c(0,0),b4=c(-1,0),summerOnly=FALSE,dayOnly=FALSE,ns,sVals){
  print(ns)
  n.cores <- length(sVals)
  registerDoParallel(cores=n.cores)
  
  dataDirectory <- "data/"
  siteData <- read.csv('/projectnb/dietzelab/kiwheel/chlorophyllCycling/allPhenocamDBsitesComplete.csv',header=TRUE)
  
  sites <- c("harvard","umichbiological","bostoncommon","coweeta","howland2",
             "morganmonroe","missouriozarks","queens","dukehw","lacclair","bbc1","NEON.D08.DELA.DP1.00033",
             "arbutuslake","bartlettir","proctor","oakridge1","hubbardbrook","canadaOA","alligatorriver","readingma",
             "bullshoals","thompsonfarm2N","ashburnham","shalehillsczo")
  yearsRemoved <- c(2015,2010,2020,2015,2017,2015,2015,2010,2017,2019,2018,2017,
                    2012,2019,2019,2010,2014,2015,2017,2018,2016,2011,2012,2019)
  nchain=5
  
  foreach(s =sVals) %dopar% {
    siteName <- sites[s]
    print(siteName)
    yearRemoved <- yearsRemoved[s]
    load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData"))
    
    variableNames <- c("p.PC","x","p.proc","b0","b4")
    inits <- list()
    for(i in 1:nchain){
      inits[[i]] <- list(b0 = min(runif(1,b0[1],b0[2]),-0.001),
                         b4 = min(runif(1,b4[1],b4[2]),-0.001),
                         p.PC = rgamma(1,1.56,0.016),
                         p.proc = rgamma(1,1.56,0.016))
    }

    includeB <- "b3"
    variableNames <- c(variableNames,"b3")
    for(n in ns){
      print(n)
      inputFileName <- paste0(siteName,"_meanTemp_summer",n,"_expBreak_slope_",includeB,"_calibration_varBurn.RData")
      if(file.exists(inputFileName)){
        load(inputFileName) #Load Model Output
        out.mat <- data.frame(as.matrix(out.burn$param))
        for(i in 1:nchain){
          inits[[i]]$b0 <- min(rnorm(1,mean(out.mat$b0),sd(out.mat$b0)),-0.001)
          inits[[i]]$b3 <- max(rnorm(1,mean(out.mat$b3),sd(out.mat$b3)),0.001)
          inits[[i]]$b4 <- min(rnorm(1,mean(out.mat$b4),sd(out.mat$b4)),-0.001)
        }
      }else{
        for(i in 1:nchain){
          inits[[i]]$b3 <- runif(1,b3[1],b3[2])
        }
      }
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
    Tair[i,yr] ~ dnorm(TairMuFinal[i,yr],TairPrecFinal[i,yr])
    
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
    
  }
    "#}
      
      if(summerOnly){
        if(dayOnly){
          outputFileName <- paste0(siteName,"_dayTemp_summer",n,"_expBreak_slope_forecast_",includeB,"_calibration_varBurn.RData")
        }else{
          outputFileName <- paste0(siteName,"_meanTemp_summer",n,"_expBreak_slope_forecast_",includeB,"_calibration_varBurn.RData")
        }
        dataFinal$p[(n+1):nrow(dataFinal$p),] <- NA
      }else{
        if(dayOnly){
          outputFileName <- paste0(siteName,"_dayTemp_fall_expBreak_slope_forecast_",includeB,"_calibration_varBurn.RData")

        }else{
          outputFileName <- paste0(siteName,"_meanTemp_fall_expBreak_slope_forecast_",includeB,"_calibration_varBurn.RData")
        }
      }
      partialFileName <- paste0("partial_",outputFileName)
      if(!file.exists(partialFileName)){
        #Remove year
        yearInt <- which(dataFinal$years==yearRemoved)
        dataFinal$p[,yearInt] <- NA
        
        ##Add priors
        dataFinal$s1.PC <- 1.56
        dataFinal$s2.PC <- 0.016
        dataFinal$s1.proc <- 1.56
        dataFinal$s2.proc <- 0.016
        dataFinal$b0_lower <- b0[1]
        dataFinal$b0_upper <- b0[2]
        dataFinal$b1_lower <- b1[1]
        dataFinal$b1_upper <- b1[2]
        dataFinal$b2_lower <- b2[1]
        dataFinal$b2_upper <- b2[2]
        dataFinal$b3_lower <- b3[1]
        dataFinal$b3_upper <- b3[2]
        dataFinal$b4_lower <- b4[1]
        dataFinal$b4_upper <- b4[2]
        
        if(dayOnly){
          dataFinal$TairMuFinal <- dataFinal$TairMuDay
          dataFinal$TairPrecFinal <- dataFinal$TairPrecDay
        }else{
          dataFinal$TairMuFinal <- dataFinal$TairMu
          dataFinal$TairPrecFinal <- dataFinal$TairPrec
        }
        
        save(file=initsFileName,inits) #Need to save inits for dic calculations
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
}
b0=c(-1,0)
b1=c(0,0)
b2=c(0,0)
b3=c(0,0)
b4=c(0,0)

createChlorophyllCyclingModelCalibration(b3=c(0,1),summerOnly=TRUE,
                                         dayOnly=FALSE,
                                         ns=seq(114,105,-1),
                                         sVals=seq(1,24)[-3]) #Change for number of included days (ns) and sites (sVals)

