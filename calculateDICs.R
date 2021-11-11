#Calculate DICs
##Create Calibration Fits for mean air temperature no trigger model
library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(suncalc)
library(rnoaa)
library(doParallel)

createDICs <- function(includeB,summerOnly=FALSE,dayOnly=FALSE,n=183,b0=c(-1,0),b1=c(0,0),b2=c(0,0),b3=c(0,0),b4=c(-1,0)){
  #n.cores <- 12
  #registerDoParallel(cores=n.cores)
  
  dataDirectory <- "data/"

  sites <- c("harvard","umichbiological","bostoncommon","coweeta","howland2",
             "morganmonroe","missouriozarks","queens","dukehw","lacclair","bbc1","NEON.D08.DELA.DP1.00033",
             "arbutuslake","bartlettir","proctor","oakridge1","hubbardbrook","asa","canadaOA","alligatorriver","readingma",
             "bullshoals","thompsonfarm2N","ashburnham","shalehillsczo")
  yearsRemoved <- c(2015,2010,2020,2015,2017,2015,2015,2010,2017,2019,2018,2017,
                    2012,2019,2019,2010,2014,2012,2015,2017,2018,2016,2011,2012,2019)
  nchain=5
  #foreach(s =1:length(sites)) %dopar% {
  for(s in 1:length(sites)){
    siteName <- sites[s]
    yearRemoved <- yearsRemoved[s]
    load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData"))
    
    if(includeB=="b1"){
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
      
      xmu[i,yr] <- max(min((x[(i-1),yr] + b4 * x[(i-1),yr]) + max(0,(b0 + b1 * Tair[i,yr])),x[1,yr]),0)
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
      b4 ~ dunif(b4_lower,b4_upper)
      }
      "
    }else if(includeB=="b2"){
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
      
      xmu[i,yr] <- max(min((x[(i-1),yr] + b4 * x[(i-1),yr]) + max(0,(b0 + b2 * D[i,yr])),x[1,yr]),0)
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
      b2 ~ dunif(b2_lower,b2_upper)
      b4 ~ dunif(b4_lower,b4_upper)
      }
      "
    }else if(includeB=="b3"){
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
    x[1,yr] ~ dbeta(x1.a[yr],x1.b[yr])
    }
    p.PC ~ dgamma(s1.PC,s2.PC)
    p.proc ~ dgamma(s1.proc,s2.proc)
    
    b0 ~ dunif(b0_lower,b0_upper)
    b3 ~ dunif(b3_lower,b3_upper)
    b4 ~ dunif(b4_lower,b4_upper)
    }
    "
    }
    
    if(summerOnly){
      if(dayOnly){
        outputFileName <- paste0(siteName,"_dayTemp_summer",n,"_expBreak_slope_",includeB,"_calibration_varBurn.RData")
        initsFileName <- paste0(siteName,"_dayTemp_summer",n,"_expBreak_slope_",includeB,"_calibration_inits.RData")
        dicFileName <- paste0(siteName,"_dayTemp_summer",n,"_expBreak_slope_",includeB,"_calibration_dic.RData")
      }else{
        outputFileName <- paste0(siteName,"_meanTemp_summer",n,"_expBreak_slope_",includeB,"_calibration_varBurn.RData")
        initsFileName <- paste0(siteName,"_meanTemp_summer",n,"_expBreak_slope_",includeB,"_calibration_inits.RData")
        dicFileName <- paste0(siteName,"_meanTemp_summer",n,"_expBreak_slope_",includeB,"_calibration_dic.RData")
      }
      dataFinal$n <- n #Limit days to July 1 through September 15th 
    }else{
      if(dayOnly){
        outputFileName <- paste0(siteName,"_dayTemp_fall_expBreak_slope_",includeB,"_calibration_varBurn.RData")
        initsFileName <- paste0(siteName,"_dayTemp_fall_expBreak_slope_",includeB,"_calibration_inits.RData")
        dicFileName <- paste0(siteName,"_dayTemp_fall_expBreak_slope_",includeB,"_calibration_dic.RData")
      }else{
        outputFileName <- paste0(siteName,"_meanTemp_fall_expBreak_slope_",includeB,"_calibration_varBurn.RData")
        initsFileName <- paste0(siteName,"_meanTemp_fall_expBreak_slope_",includeB,"_calibration_inits.RData")
        dicFileName <- paste0(siteName,"_meanTemp_fall_expBreak_slope_",includeB,"_calibration_dic.RData")
      }
    }
    if(file.exists(outputFileName) & file.exists(initsFileName)){
      if(!file.exists(dicFileName)){
        load(initsFileName) #inits
        load(outputFileName) #out.burn
        
        #Remove year
        yearInt <- which(dataFinal$years==yearRemoved)
        dataFinal$p <- dataFinal$p[,-yearInt]
        dataFinal$TairMu <- dataFinal$TairMu[,-yearInt]
        dataFinal$TairPrec <- dataFinal$TairPrec[,-yearInt]
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
        
        j.model   <- jags.model(file = textConnection(generalModel),
                                data = dataFinal,
                                inits = inits,
                                n.chains = nchain,
                                n.adapt = 1500)
        
        var.sum <- summary(out.burn$params)
        DIC <- dic.samples(j.model,n.iter = var.sum$end)
        save(DIC,file=dicFileName)
      }
    }
  }
}
createDICs(includeB="b1",summerOnly=FALSE,dayOnly=FALSE,n=183,b1=c(0,1))
createDICs(includeB="b3",summerOnly=FALSE,dayOnly=FALSE,n=183,b3=c(0,1))
  
  
  