library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(doParallel)
library('ecoforecastR')
source('generalVariables.R')

dataDirectory <- "data/"

allSites <- as.character(siteData$siteName)

allTrans <- matrix(nrow=length(allSites),ncol=14)
colnames(allTrans) <- c("siteName","meanDOY",seq(1,12))
for(s in 1:length(allSites)){
  siteName <- allSites[s]
  print(siteName)
  transDates <- numeric()
  load(paste0(dataDirectory,siteName,"_dataFinal.RData"))
  for(yr in 1:dataFinal$N){
    yrName <- dataFinal$years[yr]
    print(yrName)
    p.file <- paste0(transitionEstimateOutputsFolder,siteName,"_",yrName,"_PhenoCam_changePointCurve_varBurn.RData")
    if(file.exists(p.file)){
      load(p.file)
      if(typeof(var.burn)!=typeof(FALSE)){
        var.mat <- data.frame(as.matrix(var.burn))

        if(sd(var.mat$k)<7){
          transDates <- c(transDates,(182+mean(var.mat$k)))
        }else{
          transDates <- c(transDates,NA)
        }
      }else{
        transDates <- c(transDates,NA)
      }
    }else{
      transDates <- c(transDates,NA)
    }
  }
  allTrans[s,2] <- round(mean(na.omit(transDates)),digits=2)
  allTrans[s,3:(dataFinal$N+2)] <- round(transDates - as.numeric(allTrans[s,2]),digits=2)
  allTrans[s,1] <- siteName
}
write.table(x=allTrans,sep=",",file=allPhenoTranFile,row.names=FALSE,quote = FALSE)
