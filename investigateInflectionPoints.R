library('rjags')
library('runjags')
library('RColorBrewer')
library(doParallel)
source('generalVariables.R')
source('calculateStart.R')

##Investigate model predictions!
registerDoParallel(cores=n.cores)
ns <- seq(35,183)

transData <- read.csv(allPhenoTranFile,header=TRUE)

foreach(s =1:length(sites)) %dopar% {
  outputData <- matrix(nrow=0,ncol=8)
  colnames(outputData) <- c("siteName","year","includedYear","trans","n","converged","minDiffDiff","meanDiffDiff")
  siteName <- as.character(sites[s])
  print(siteName)

  load(paste0(dataDirectory,siteName,"_dataFinal.RData"))
  siteTrans <- transData[transData$siteName==siteName,]
  for(n in ns){
    fileName <- paste0(CCmodelOutputsFolder,siteName,"_",n,"_ccModel_forecast_calibration_varBurn.RData")
    if(file.exists(fileName)){
      res <- try(load(fileName))#Load Model Output 
      if(inherits(res,"try-error")){
        next
      }
    }else{
      next 
    }
    pred.mat <- data.frame(as.matrix(out.burn$predict))
    
    prow = sample.int(nrow(pred.mat),Nmc,replace=TRUE)
    pred.mat <- pred.mat[prow,]
    out.mat <- data.frame(as.matrix(out.burn$param))
    
    b0 <- out.mat$b0[prow]
    b4 <- out.mat$b4[prow]
    b3 <- out.mat$b3[prow]
    if(diff(range(b0))<0.5 & diff(range(b4))<0.5 & diff(range(b3))<0.5){
      converged <- TRUE
    }else{
      converged <- FALSE
    }

    for(yr in 1:dataFinal$N){
      yrName <- dataFinal$years[yr]
      print(yrName)
      if(converged){
        yrPred <- pred.mat[1:Nmc,(184*(yr-1)+1):(184*(yr-1)+184)]
        
        avgDiffs <- matrix(nrow=Nmc,ncol=ncol(yrPred))
        for(g in 1:Nmc){
          avgDiffs[g,] <- calculateStart(ys=as.numeric(yrPred[g,]),avgNum = avgNum)
        }
        ci <- apply(avgDiffs,quantile,MARGIN = 2,c(0.025,0.5,0.975))
        pred.ci <- apply(yrPred,quantile,MARGIN = 2,c(0.025,0.5,0.975))
        newRow <- c(siteName,yrName,yrName!=dataFinal$yearRemoved,(siteTrans[yr+2]+siteTrans[2]),
                    n,converged,round(min(ci),digits=4),round(min(ci[2,]),digits=4))
      }else{
        newRow <- c(siteName,yrName,yrName!=dataFinal$yearRemoved,(siteTrans[yr+2]+siteTrans[2]),n,converged,NA,NA)
      }
      outputData <- rbind(outputData,newRow)
    }
  }
  outputData <- data.frame(outputData)
  outputData$offset <- round((as.numeric(outputData$n)+182)-as.numeric(outputData$trans),digits=0)
  
  save(outputData,file=paste0(siteName,"_inflectionPointData_15.RData"))
}

