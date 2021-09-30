library('rjags')
library('runjags')
library('ecoforecastR')
library('PhenologyBayesModeling')

#Create forecast step model:
#' Basic logistic forecast step
#'
#' @param IC Initial conditions
#' @param b0 The parameter b0
#' @param b1 The parameter b1
#' @param b2 The parameter b2
#' @param b3 The parameter b3
#' @param Q Process error (default = 0 for deterministic runs)
#' @param n Size of Monte Carlo ensemble
#' @param NT number of days for forecast
#' @param D
#' @import rjags
#' @import runjags
#' @import ecoforecastR
#'
#' @return
#' @export
#'
#' @examples
forecastStep <- function(IC,b0,b1,b2,b3,Q=0,n,NT,Tair,D){
  x <- matrix(NA,n,NT)
  if(length(IC)==1){
    Xprev <- rep(IC,n)
  }else{
    Xprev <- IC
  }
  
  for(t in 1:NT){
    bd <- Xprev + b0 + (b1 * Xprev) + (b2 * Xprev **2) 
    syn <- b3 * Tair[t]
    if(length(syn)==1){
      syn <- rep(syn,n)
    }
    
    xNew <- numeric()
    mu <- rep(NA,n)
    for(i in 1:n){
      #mu[i] <- max(0,min(mu[i],0.999))
      mu[i] <- bd[i] + max(syn[i],0)
      mu[i] <- max(0,min(mu[i],IC[min(i,length(IC))]))
      
      if(length(Q)>1){
        xNew <- c(xNew,rnorm(1,mu[i],Q[i]))
      }else{
        xNew <- c(xNew,rnorm(1,mu[i],Q))
      }
      #xNew <- c(xNew,max(0, min(0.99,xl)))
    }
    x[,t] <- xNew ## trunate normal process error
    Xprev <- x[,t]
  }
  return(x)
}

calTransRoots <- function(c,b,m,d){
  root1 <- (b*m+log(sqrt(3)+2))/b
  root2 <- (b*m+log(2-sqrt(3)))/b
  return(c(root1,root2))
}

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
load('harvard_b3_final_calibration_varBurn.RData')

Nmc <- 1000
out.mat <- data.frame(as.matrix(out.burn$param))

b0 <- out.mat$b0
b1 <- out.mat$b1
b2 <- out.mat$b2
b3 <- out.mat$b3
prow = sample.int(nrow(out.mat),Nmc,replace=TRUE)
baseTempOrig <- 20
sites <- c("asuhighlands","bullshoals","macleish","harvard")

#files <- dir(pattern = "phenoCurve_varBurn.RData")
files <- dir(pattern = "_ysDet_changePointCurve_varBurn.RData")
pdf(file="PhenoFits_transitionEstimates.pdf",height=6,width=10)

for (f in 1:length(files)){
  siteName <- strsplit(files[f],"_")[[1]][1]
  print(siteName)
  yr <- as.numeric(strsplit(files[f],"_")[[1]][2])
  
  load(files[f])
  
  fileName <- paste0(siteName,"_b3_final_calibration_varBurn.RData")
  
  load(fileName)
  load(paste(siteName,"_CDD",baseTempOrig,"_dataFinal.RData",sep=""))
  
  
  out.mat <- data.frame(as.matrix(out.burn$param))
  NT <- length(dataFinal$TairMu[,1])
  out.mat.pred <- data.frame(as.matrix(out.burn$predict))
  
  prow.pred = sample.int(nrow(out.mat),Nmc,replace=TRUE)
  initialXs <- out.mat.pred[prow.pred,1]
  days <- seq(1,NT)
  
  initialXs <- out.mat.pred[prow.pred,dataFinal$n*(yr-1)+1]
  ysDet <- as.numeric(forecastStep(IC=mean(initialXs),b0=mean(b0),b1=mean(b1),b2=mean(b2),
                                   b3=mean(b3),n=1,NT=length(days),Tair=dataFinal$TairMu[,yr],
                                   D=dataFinal$D[,yr]))
  
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
  # TranF <- var.mat$TranF[rndNums]
  # bF <- var.mat$bF[rndNums]
  # c <- var.mat$c[rndNums]
  # d <- rep(0,length(rndNums))
  # for(g in 1:10000){
  #   ycred[g,] <- pheno.logistic(Tran=TranF[g],b=bF[g],c=c[g],d=d[g],xseq=days)
  #   tranCred1[g] <- calTransRoots(c=c[g],b=bF[g],m=TranF[g],d=d[g])[1]
  #   tranCred2[g] <- calTransRoots(c=c[g],b=bF[g],m=TranF[g],d=d[g])[2]
  # }
  ci <- apply(ycred,2,quantile,c(0.025,0.5,0.975))
  tran.ci <- quantile(k,c(0.025,0.5,0.975))
  # tran1.ci <- quantile(tranCred1,c(0.025,0.5,0.975))
  # tran2.ci <- quantile(tranCred2,c(0.025,0.5,0.975))
  
  ##PhenoCam
  p.file <- paste0(siteName,"_",yr,"_p_phenoCurve_varBurn00000.RData")
  if(file.exists(p.file)){
    load(p.file)
    p.var.burn <- var.burn
    var.mat <- data.frame(as.matrix(p.var.burn))
    ycred <- matrix(0,nrow=10000,ncol=length(days))
    tranCred1 <- rep(0,10000)
    tranCred2 <- rep(0,10000)
    
    rndNums <- sample(1:nrow(var.mat),10000,replace=T)
    TranF <- var.mat$TranF[rndNums]
    bF <- var.mat$bF[rndNums]
    c <- var.mat$c[rndNums]
    d <- rep(0,length(rndNums))
    for(g in 1:10000){
      ycred[g,] <- pheno.logistic(Tran=TranF[g],b=bF[g],c=c[g],d=d[g],xseq=days)
      tranCred1[g] <- calTransRoots(c=c[g],b=bF[g],m=TranF[g],d=d[g])[1]
      tranCred2[g] <- calTransRoots(c=c[g],b=bF[g],m=TranF[g],d=d[g])[2]
    }
    p.ci <- apply(ycred,2,quantile,c(0.025,0.5,0.975))
    p.tran1.ci <- quantile(tranCred1,c(0.025,0.5,0.975))
    p.tran2.ci <- quantile(tranCred2,c(0.025,0.5,0.975))
  }
  
  plot(days,ysDet,type="l",main=paste(siteName,yr))
  
  ecoforecastR::ciEnvelope(days,ci[1,],ci[3,],col=col.alpha("lightblue",0.5))
  if(file.exists(p.file)){
    ecoforecastR::ciEnvelope(days,p.ci[1,],p.ci[3,],col=col.alpha("pink",0.5))
  }
  lines(days,ysDet)
  points(days,dataFinal$p[,yr],pch=20)
  #abline(v=tran1.ci[2],col="red")
  #abline(v=tran1.ci[1],col="red",lty=2)
  #abline(v=tran1.ci[3],col="red",lty=2)
  
  abline(v=tran.ci[2],col="blue")
  abline(v=tran.ci[1],col="blue",lty=2)
  abline(v=tran.ci[3],col="blue",lty=2)
  if(file.exists(p.file)){
    abline(v=p.tran2.ci[2],col="red")
    abline(v=p.tran2.ci[1],col="red",lty=2)
    abline(v=p.tran2.ci[3],col="red",lty=2)
  }
}
dev.off()















