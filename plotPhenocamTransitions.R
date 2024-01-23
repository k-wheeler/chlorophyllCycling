library(rjags)
library(runjags)
source('ciEnvelope.R')
source('generalVariables.R')

changepointModel <- function(y1,mS,mF,k,xseq){
  yseq <- y1
  for(x in xseq){
    if(x<=k){
      yseq <- c(yseq,yseq[length(yseq)]-mS)
    }else{
      yseq <- c(yseq,yseq[length(yseq)]-mF)
    }
  }
  return(yseq)
}

pdf(file="PhenoCam_transitionEstimates_updatedK.pdf",height=6,width=10)

for(s in 1:nrow(siteData)){
  siteName <- as.character(siteData$siteName[s])
  print(siteName)
  if(file.exists(paste0(dataDirectory,siteName,"_dataFinal.RData"))){
    load(paste0(dataDirectory,siteName,"_dataFinal.RData"))
    days <- seq(1,dataFinal$n)
    for(yr in 1:dataFinal$N){
      if(length(which(dataFinal$p[,yr]<0.10))>0){
        p <- dataFinal$p[1:(which(dataFinal$p[,yr]<0.10)[1]-1),yr]
      }else{
        p <- dataFinal$p[,yr]
      }
      yrName <- dataFinal$years[yr]
      print(yrName)
      p.file <- paste0(transitionEstimateOutputsFolder,siteName,"_",yrName,"_PhenoCam_changePointCurve_varBurn.RData")
      plot(days,dataFinal$p[,yr],pch=20,main=paste(siteName,yrName))
      abline(h=0.15,col="red")
      
      if(file.exists(p.file)){
        load(p.file)
        if(typeof(var.burn)!=typeof(FALSE)){
          var.mat <- data.frame(as.matrix(var.burn))
          ycred <- matrix(0,nrow=10000,ncol=length(days))
          
          rndNums <- sample(1:nrow(var.mat),10000,replace=T)
          mS <- var.mat$mS[rndNums]
          mF <- var.mat$mF[rndNums]
          y1 <- var.mat$y.1[rndNums]
          k <- var.mat$k[rndNums]
          for(g in 1:10000){
            ycred[g,] <- changepointModel(y1=y1[g],mS=mS[g],mF=mF[g],k=k[g],xseq=days[2:length(days)])
          }
          ci <- apply(ycred,2,quantile,c(0.025,0.5,0.975))
          tran.ci <- quantile(k,c(0.025,0.5,0.975))
          
          ciEnvelope(days,ci[1,],ci[3,],col=col.alpha("lightblue",0.5))
          points(days,dataFinal$p[,yr],pch=20)
          abline(v=tran.ci[2],col="blue")
          abline(v=tran.ci[1],col="blue",lty=2)
          abline(v=tran.ci[3],col="blue",lty=2)
        }
      }
    }
  }
}
dev.off()
