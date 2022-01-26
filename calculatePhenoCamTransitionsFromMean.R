library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(doParallel)
library('ecoforecastR')

dataDirectory <- "data/"
sites <- c("harvard","umichbiological","bostoncommon","coweeta","howland2",
           "morganmonroe","missouriozarks","queens","dukehw","lacclair","bbc1","NEON.D08.DELA.DP1.00033",
           "arbutuslake","bartlettir","proctor","oakridge1","hubbardbrook","asa","canadaOA","alligatorriver","readingma",
           "bullshoals","thompsonfarm2N","ashburnham","shalehillsczo")

siteData <- read.csv('/projectnb/dietzelab/kiwheel/chlorophyllCycling/allPhenocamDBsitesComplete.csv',header=TRUE)
sites <- as.character(siteData$siteName)

allTrans <- matrix(nrow=length(sites),ncol=14)
colnames(allTrans) <- c("siteName","meanDOY",seq(1,12))
for(s in 1:length(sites)){
  siteName <- sites[s]
  print(siteName)
  transDates <- numeric()
  load(paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData"))
  for(yr in 1:dataFinal$N){
    yrName <- dataFinal$years[yr]
    print(yrName)
    p.file <- paste0(siteName,"_",yrName,"_PhenoCam_changePointCurve_varBurn.RData")
    if(file.exists(p.file)){
      load(p.file)
      if(typeof(var.burn)!=typeof(FALSE)){
        var.mat <- data.frame(as.matrix(var.burn))
        ycred <- matrix(0,nrow=10000,ncol=length(days))
        
        transDates <- c(transDates,(182+mean(var.mat$k)))
        
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
write.table(x=allTrans,sep=",",file="phenocamTransitions_fromMean.csv",row.names=FALSE,quote = FALSE)
