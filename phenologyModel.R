#' General phenology forecast function
#'
#' @param siteName The site name
#' @param URLs The URL where the site's phenocam data is located
#' @param lat The latitude of the site in decimals
#' @param long The longitude of the site in decimals
#' @param dataDirectory The file path for the directory to download and store data in
#' @param startDate The start date for the forecast in date form
#' @param endDate The end date for the forecast in date form
#' @param baseTemp If CDD, the desired base temperature
#' @param ERA5dataFolder
#' @param baseNum
#' @param iterSize
#' @param effSize
#' @param partialFile
#' @param vars
#' @param TZ
#' @param fittedDat
#' @param splitYears
#' @import rjags
#' @import runjags
#' @import suncalc
#' @import phenopix
#' @import coda
#' @import PhenologyBayesModeling
#' @export
phenologyCalibration_Autumn <- function(siteName,
                                        URLs,lat,long,dataDirectory,startDate,endDate,
                                        baseTemp=NA,ERA5dataFolder,
                                        baseNum=10000,
                                        iterSize=5000,
                                        effSize=5000,
                                        vars,TZ,calDIC=FALSE,
                                        fittedDat, splitYears=""){

nchain=5
###Download PhenoCam data and format
phenoData <- matrix(nrow=0,ncol=32)
print(URLs[1])
for(u in 1:length(URLs)){
  phenoDataSub <- download.phenocam(URLs[u])
  phenoData <- rbind(phenoData,phenoDataSub)
}
#print(phenoData)
##Order and remove duplicate PC data
phenoData2 <- phenoData[order(phenoData$date),]
phenoData3 <- phenoData2[!duplicated(phenoData2$date),]
phenoData <- phenoData3

phenoData <- phenoData[phenoData$date<endDate,]
p.old <- phenoData$gcc_90
time.old <-  as.Date(phenoData$date)
#print("passed time.old")

days <- seq(as.Date(startDate),(as.Date(endDate)),"day")
#print("past days")
p <- rep(NA,length(days))

for(i in 1:length(p.old)){
  p[which(days==time.old[i])] <- p.old[i]
}

months <- lubridate::month(days)
years <- lubridate::year(days)

dat2 <- data.frame(dates=days,years=years,months=months,p=p)
datTairEns <- load_ERA5_Tair_New(ERA5dataFolder=ERA5dataFolder,endDate=endDate)
print("Finished loading ERA5")
TairMu <- apply(X=datTairEns,MARGIN=2,FUN=mean)
TairPrec <- 1/apply(X=datTairEns,MARGIN=2,FUN=var)
dat2$TairMu <- TairMu # + (0-min(TairMu)) ##Done to make sure all temperatures are >= 0 
#baseTempOrig <- baseTemp
#baseTemp <- baseTemp + (0-min(TairMu)) 
dat2$TairPrec<- TairPrec

if(TZ==5){
  TZ_name <- "America/New_York"
}else if (TZ==6){
  TZ_name <- "America/Chicago"
}else{
  print("Unrecognized TZ")
}

dayLengths <- numeric()

for(d in 1:length(days)){
  suntimes <- getSunlightTimes(date=days[d],
                               lat=lat,lon=long,keep=c("nauticalDawn","nauticalDusk"),
                               tz = TZ_name)
  dayLengths <- c(dayLengths,as.numeric(suntimes$nauticalDusk-suntimes$nauticalDawn))
}

dat2$D <- dayLengths
ICsdat <- dat2[as.numeric(format(dat2$dates,"%j"))%in% seq(203,212),]
dat2 <- dat2[as.numeric(format(dat2$dates,"%j"))%in% seq(213,365),]

nrowNum <- 365-212
p <- matrix(nrow=nrowNum,ncol=0)
TairMu <- matrix(nrow=nrowNum,ncol=0)
D <- matrix(nrow=nrowNum,ncol=0)
ICs <- matrix(nrow=10,ncol=0)
TairPrec <- matrix(nrow=nrowNum,ncol=0)
valNum <- 0
days2 <- matrix(nrow=nrowNum,ncol=0)

finalYrs <- numeric()
sofs <- numeric()
for(i in (lubridate::year(as.Date(dat2$dates[1]))):lubridate::year(as.Date(dat2$dates[length(dat2$dates)]))){##I know this includes the forecasted stuff, but it shouldn't really matter because of the JAGS model setup
  #print(i)
  subDat <- dat2[lubridate::year(as.Date(dat2$dates))==i,]
  valNum <- valNum + 1
  Low <- fittedDat[valNum,'Low']
  High <- fittedDat[valNum,'High']
  if(!is.na(Low)){
    newCol <- scales::rescale(subDat$p,to=c(0,1),from=c(Low,High))
    p <- cbind(p,newCol)
    ICs <- cbind(ICs,scales::rescale(ICsdat[lubridate::year(as.Date(ICsdat$dates))==i,]$p,from=c(Low,High)))
    days2 <- cbind(days2,as.Date(subDat$dates))
    finalYrs <- c(finalYrs,i)
    sofs <- c(sofs,(fittedDat[valNum,'FallStartDay']-212)) ######Change for start if needed
    TairMu <- cbind(TairMu,subDat$TairMu)
    D <- cbind(D,subDat$D)
    TairPrec <- cbind(TairPrec,subDat$TairPrec)
  }
}
p[p<0] <- 0
p[p>0.999] <- 0.999
ICs[ICs<0] <- 0
ICs[ICs>0.999] <- 0.999

rndNums <- seq(1,floor(ncol(p)/2))
rndNums2 <- seq(ncol(p)-length(rndNums)+1,ncol(p))
print(paste("Initial dim(p):",dim(p)))
if(splitYears=="first"){
  p <- p[,rndNums]
  TairMu <- TairMu[,rndNums]
  TairPrec <- TairPrec[,rndNums]
}else if(splitYears=="second"){
  p <- p[,rndNums2]
  TairMu <- TairMu[,rndNums2]
  TairPrec2 <- TairPrec[,rndNums2]
}
print(paste("Final dim(p):",dim(p)))

dataFinal <- list(p=p,years=finalYrs,sofMean=mean(sofs))
dataFinal$n <- nrowNum
dataFinal$N <- ncol(dataFinal$p)
dataFinal$CDDtrigger.lower <- 0
dataFinal$CDDtrigger.upper <- 500
dataFinal$s1.PC <- 1.56#1262.626 ## Very roughly based off of what I think are reasonable and uninformed priors
dataFinal$s2.PC <- 0.016#50.50505 
dataFinal$s1.proc <- 1.56#1262.626
dataFinal$s2.proc <- 0.016#50.50505

x1a <- numeric()
x1b <- numeric()
for(yr in 1:dataFinal$N){
  mu <- mean(ICs[,yr],na.rm=TRUE)
  vr <- var(ICs[,yr],na.rm = TRUE)
  x1a <- c(x1a,(mu**2-mu**3-mu*vr)/(vr))
  x1b <- c(x1b,(mu-2*mu**2+mu**3-vr+mu*vr)/(vr))
}

dataFinal$x1.a <- x1a
dataFinal$x1.b <- x1b

dataFinal$TairMu <- TairMu
dataFinal$TairPrec <- TairPrec
dataFinal$baseTemp <- baseTemp
dataFinal$D <- D

inits <- list()
save(dataFinal,file=paste(siteName,"_CDD",baseTemp,"_dataFinal.RData",sep=""))
if(vars=="noCov"){
  outputFileName <- paste(siteName,"_CDD20_final_calibration_varBurn.RData",sep="")
  
  variableNames <- c("p.PC","p.proc","x","b0","b1","b2","a","CDDtrigger")
  dataFinal$b0_lower <- -1
  dataFinal$b0_upper <- 1#1#0
  dataFinal$b1_lower <- -1#-1
  dataFinal$b1_upper <- 1
  dataFinal$b2_lower <- -1
  dataFinal$b2_upper <- 1#1
  dataFinal$a_upper <- 0
  dataFinal$a_lower <- -0.01
  dataFinal$CDDtrigger.lower <- 0
  dataFinal$CDDtrigger.upper <- 500
  
  generalModel = "
    model {
  ### Data Models
  for(yr in 1:(N)){
  for(i in 1:n){
  p[i,yr] ~ dnorm(x[i,yr],p.PC)
  }
  }
  
  #### Process Model
  for(yr in 1:(N)){
  for(i in 2:n){
  Tair[i,yr] ~ dnorm(TairMu[i,yr],TairPrec[i,yr])
  CDDs[i,yr] <- ifelse(Tair[i,yr]<baseTemp,CDDs[(i-1),yr]+baseTemp - Tair[i,yr],CDDs[(i-1),yr])
  xmu[i,yr] <- max(min(x[(i-1),yr] + ifelse(CDDs[i,yr]>CDDtrigger,(b0 + (b1 * x[(i-1),yr]) + (b2 * x[(i-1),yr] ** 2)),a),x[1,yr]),0)
  x[i,yr] ~ dnorm(xmu[i,yr],p.proc) # T(0,0.999)
  }
  }
  
  #### Priors
  for(yr in 1:N){ ##Initial Conditions
  x[1,yr] ~ dbeta(x1.a[yr],x1.b[yr])
  CDDs[1,yr] <- 0
  }
  p.PC ~ dgamma(s1.PC,s2.PC)
  p.proc ~ dgamma(s1.proc,s2.proc)
  CDDtrigger ~ dunif(CDDtrigger.lower,CDDtrigger.upper)
  
  a ~ dunif(a_lower,a_upper)
  b0 ~ dunif(b0_lower,b0_upper)
  b2 ~ dunif(b2_lower,b2_upper)
  b1 ~ dunif(b1_lower,b1_upper) 
}
"

  ##Determine Inits:
  transitionCDD <- numeric()
  agings <- numeric()
  for(yr in 1:dataFinal$N){
    CDD <- 0
    for(i in 1:sofs[yr]){
      if(dataFinal$TairMu[i,yr]<dataFinal$baseTemp){
        CDD <- CDD + (dataFinal$baseTemp - dataFinal$TairMu[i,yr])
      }
    }
    transitionCDD <- c(transitionCDD,CDD)
    agings <- c(agings,
                as.numeric(lm(dataFinal$p[,yr]~seq(1,length(dataFinal$p[,yr])))$coefficients[2]))
  }
  
  for(i in 1:nchain){
    inits[[i]] <- list(CDDtrigger = rnorm(1,mean(transitionCDD),10),
                       a=rnorm(1,mean(agings),0.0005),
                       b0=rnorm(1,-0.25,0.05),
                       b1=rnorm(1,0.6,0.08),
                       b2=rnorm(1,-0.5,0.1))
  }
}else if(vars=="b3"){
  outputFileName <- paste0(siteName,"_b3_final_calibration_varBurn.RData")
  dataFinal$b0_lower <- -1
  dataFinal$b0_upper <- 0
  dataFinal$b1_lower <- -1
  dataFinal$b1_upper <- 0
  dataFinal$b2_lower <- -1
  dataFinal$b2_upper <- 0
  dataFinal$b3_lower <- 0
  dataFinal$b3_upper <- 1
  
  variableNames <- c("p.PC","x","b0","b1","b2","b3","p.proc")
  
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
  
  xmu[i,yr] <- max(min(x[(i-1),yr] + (b0 + (b1 * x[(i-1),yr]) + (b2 * x[(i-1),yr] ** 2)) + max(0,(b3 * TairMu[i,yr])),x[1,yr]),0)
  x[i,yr] ~ dnorm(xmu[i,yr],p.proc) #T(0,0.999)
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
  b3 ~ dunif(b3_lower,b3_upper)
  }
  "
}

j.model   <- jags.model(file = textConnection(generalModel),
                        data = dataFinal,
                        inits = inits,
                        n.chains = nchain,
                        n.adapt = 1500)

out.burn <- runForecastIter(j.model=j.model,variableNames=variableNames,
                            baseNum = 15000,iterSize = 5000,effSize = 5000,partialFile = paste("partial_",outputFileName,sep=""))
##Thin the data:
out.mat <- as.matrix(out.burn$params)
thinAmount <- round(nrow(out.mat)/5000,digits=0)
out.burn2 <- list()
out.burn2$params <- window(out.burn$params,thin=thinAmount)
out.burn2$predict <- window(out.burn$predict,thin=thinAmount)
out.burn <- out.burn2
save(out.burn,file = outputFileName)
}
