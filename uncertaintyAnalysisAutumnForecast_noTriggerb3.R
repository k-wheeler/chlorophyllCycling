library('rjags')
library('runjags')
library('ecoforecastR')
library('RColorBrewer')

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

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
#' @import rjags
#' @import runjags
#' @import ecoforecastR
#'
#' @return
#' @export
#'
#' @examples
forecastStep <- function(IC,b0,b1,b2,b3,Q=0,n,NT,Tair){
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

files <- dir(pattern = "_b3_final_calibration_varBurn.RData")

#sites <- c('bostoncommon','caryinstitute','harvard','hubbardbrooksfws','proctor','marcell','umichbiological')
sites <- character()
#cls <- c('black','blue','red','green','cyan','pink','orange')
pdf(file="UncertaintyAnalysis_noTrigger_b3_calibration.pdf",height=6,width=40)

allb0s <- matrix(ncol=0,nrow=4000)
allb1s <- matrix(ncol=0,nrow=4000)
allb2s <- matrix(ncol=0,nrow=4000)
allb3s <- matrix(ncol=0,nrow=4000)

for(f in 1:length(files)){
#for(f in 1:2){
  
  fileName <- files[f]
  if(strsplit(fileName,"_")[[1]][1]!="partial"){
    print(fileName)
    tle <- strsplit(fileName,"_")[[1]][1:2]
    siteName <- strsplit(fileName,"_")[[1]][1]
    sites <- c(sites,siteName)
    # siteName <- substr(siteName,1,(nchar(siteName)-5))
    print(siteName)
    if(TRUE){
      #if(siteName %in% sites){
      #baseTempOrig <- strsplit(fileName,"_")[[1]][2]
      #baseTempOrig <- as.numeric(substr(baseTempOrig,4,5))
      baseTempOrig <- 20
      load(fileName)
      load(paste(siteName,"_CDD",baseTempOrig,"_dataFinal.RData",sep=""))
      if(ncol(dataFinal$TairMu)>1){
        out.mat <- data.frame(as.matrix(out.burn$param))
        
        b0 <- out.mat$b0
        b1 <- out.mat$b1
        b2 <- out.mat$b2
        b3 <- out.mat$b3
        NT <- length(dataFinal$TairMu[,1])
        
        out.mat.pred <- data.frame(as.matrix(out.burn$predict))
        
        Nmc <- 1000
        prow = sample.int(nrow(out.mat),Nmc,replace=TRUE)
        
        initialXs <- out.mat.pred[prow,1]
        print(mean(initialXs))
        
        days <- seq(1,NT)
        ysDet <- forecastStep(IC=mean(initialXs),b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),n=1,NT=length(days),Tair=dataFinal$TairMu[,1]) #Dependent on IC, but independent for >=1
        ysParam <- forecastStep(IC=mean(initialXs),b0=b0[prow],b1=b1[prow],b2=b2[prow],b3=b3[prow],n=Nmc,NT=NT,Tair=dataFinal$TairMu[,1])
        ysProc <- forecastStep(IC=mean(initialXs),b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),n=Nmc,NT=NT,Q=sqrt(1/out.mat$p.proc[prow]),Tair=dataFinal$TairMu[,1])
        ysIC <- forecastStep(IC=initialXs,b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),n=length(initialXs),NT=NT,Tair=dataFinal$TairMu[,1])
        
        for(i in 2:dataFinal$N){
          print(colnames(out.mat.pred)[dataFinal$n*(i-1)+1])
          initialXs <- out.mat.pred[prow,dataFinal$n*(i-1)+1]
          
          ysDet <- cbind(ysDet,forecastStep(IC=mean(initialXs),b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),n=1,NT=length(days),Tair=dataFinal$TairMu[,i]))
          ysParam <- cbind(ysParam,forecastStep(IC=mean(initialXs),b0=b0[prow],b1=b1[prow],b2=b2[prow],b3=b3[prow],n=Nmc,NT=NT,Tair=dataFinal$TairMu[,i]))
          ysProc <- cbind(ysProc,forecastStep(IC=mean(initialXs),b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),n=Nmc,NT=NT,Q=sqrt(1/out.mat$p.proc[prow]),Tair=dataFinal$TairMu[,i]))
          ysIC <- cbind(ysIC,forecastStep(IC=initialXs,b0=mean(b0),b1=mean(b1),b2=mean(b2),b3=mean(b3),n=length(initialXs),NT=NT,Tair=dataFinal$TairMu[,i]))
        }
        N.IP.ci = apply(ysParam,2,quantile,c(0.025,0.5,0.975))
        N.Proc.ci = apply(ysProc,2,quantile,c(0.025,0.5,0.975))
        N.IC.ci = apply(ysIC,2,quantile,c(0.025,0.5,0.975))
        
        # plot(days,ysDet,pch=20,ylim=c(-0.3,1.3),main=tle)
        # ecoforecastR::ciEnvelope(days,N.Proc.ci[1,],N.Proc.ci[3,],col=col.alpha("green",0.5))
        # ecoforecastR::ciEnvelope(days,N.IP.ci[1,],N.IP.ci[3,],col=col.alpha("blue",0.5))
        # ecoforecastR::ciEnvelope(days,N.IC.ci[1,],N.IC.ci[3,],col=col.alpha("pink",0.5))
        # lines(days,N.Proc.ci[2,],col="green",lwd=2)
        # lines(days,N.IP.ci[2,],col="blue",lwd=2)
        # lines(days,N.IC.ci[2,],col="pink",lwd=2)
        # lines(days,ysDetnob3,col="cyan",lwd=3)
        # legend("bottomright",c('Process','Parameter','IC'),col=c("green","blue","pink"),lty = c(1,1,1))
        # print("done first plot")
        
        # plot(density(out.mat$CDDtrigger),main=paste(tle,"CDDtrigger"))#,quantile(out.mat$CDDtrigger,c(0.25,0.5,0.75))))
        
        plot(seq(1,length(dataFinal$p)),dataFinal$p,pch=20,main=paste(tle, "p"),ylim=c(0,1))
        # N.Proc.ci <- rbind(rep(N.Proc.ci[1,],dataFinal$N),rep(N.Proc.ci[2,],dataFinal$N),rep(N.Proc.ci[3,],dataFinal$N))
        # N.IP.ci <- rbind(rep(N.IP.ci[1,],dataFinal$N),rep(N.IP.ci[2,],dataFinal$N),rep(N.IP.ci[3,],dataFinal$N))
        # N.IC.ci <- rbind(rep(N.IC.ci[1,],dataFinal$N),rep(N.IC.ci[2,],dataFinal$N),rep(N.IC.ci[3,],dataFinal$N))
        
        ecoforecastR::ciEnvelope(seq(1,length(dataFinal$p)),N.Proc.ci[1,],N.Proc.ci[3,],col=col.alpha("green",0.5))
        ecoforecastR::ciEnvelope(seq(1,length(dataFinal$p)),N.IP.ci[1,],N.IP.ci[3,],col=col.alpha("blue",0.5))
        ecoforecastR::ciEnvelope(seq(1,length(dataFinal$p)),N.IC.ci[1,],N.IC.ci[3,],col=col.alpha("pink",0.5))
        lines(seq(1,length(dataFinal$p)),N.Proc.ci[2,],col="green",lwd=2)
        lines(seq(1,length(dataFinal$p)),N.IP.ci[2,],col="blue",lwd=2)
        lines(seq(1,length(dataFinal$p)),N.IC.ci[2,],col="pink",lwd=2)
        lines(seq(1,length(dataFinal$p)),ysDet,col="cyan",lwd=3)
        #lines(days,ysDet,col="cyan",lwd=3)
        
        # pred.ci <- apply(X=out.mat.pred,2,FUN = quantile,c(0.025,0.5,0.975))
        # ecoforecastR::ciEnvelope(seq(1,ncol(pred.ci)),pred.ci[1,],pred.ci[3,],col=col.alpha("lightblue",0.5))
        # lines(seq(1,ncol(out.mat.pred)),pred.ci[2,])
        print("done second plot")
        
        #triggeredDays <- seq(min(which(dataFinal$D[,1]<13.5)),length(dataFinal$p),nrow(dataFinal$p))
        #abline(v=triggeredDays,col="red",lty=2)
        
        par(mfrow=c(1,4))
        plot(density(b0),main=paste(siteName,"b0"))
        plot(density(b1),main=paste(siteName,"b1"))
        plot(density(b2),main=paste(siteName,"b2"))
        plot(density(b3),main=paste(siteName,"b3"))
        par(mfrow=c(1,1))
        allb0s <- cbind(allb0s,b0[1:4000])
        allb1s <- cbind(allb1s,b1[1:4000])
        allb2s <- cbind(allb2s,b2[1:4000])
        allb3s <- cbind(allb3s,b3[1:4000])
      }
    }
  }
}
par(mfrow=c(1,4))
#print(head(allb0s))
cls <- seq(3,23)
plot(density(allb0s[,1]),xlim=c(-0.2,0),main="All b0",col=c25[1])
for(i in 2:ncol(allb0s)){
  lines(density(allb0s[,i]),col=c25[i])
}
plot(density(allb1s[,1]),xlim=c(-0.1,0),main="All b1")
for(i in 2:ncol(allb1s)){
  lines(density(allb1s[,i]),col=c25[i])
}
plot(density(allb2s[,1]),xlim=c(-0.1,0),main="All b2")
for(i in 2:ncol(allb2s)){
  lines(density(allb2s[,i]),col=c25[i])
}
plot(density(allb3s[,1]),xlim=c(0,1.5*max(allb3s)),main="All b3")
for(i in 2:ncol(allb3s)){
  lines(density(allb3s[,i]),col=c25[i])
}
par(mfrow=c(1,1))
plot(numeric(),numeric(),ylim=c(0,1),xlim=c(0,1))
legend("topleft",sites,col=c25,lty=rep(1,length(sites)))
dev.off()
# 
# b0Summary <- data.frame(matrix(ncol=3,nrow=ncol(allb0s)))
# colnames(b0Summary) <- c("site","mean","sd")
# b0Summary$site <- sites
# b0Summary$mean <- colMeans(allb0s)
# b0Summary$sd <- apply(X=allb0s,FUN=sd,MARGIN = 2)
# save(b0Summary,file="b0Summary_noTrigger.RData")
# 
# b1Summary <- data.frame(matrix(ncol=3,nrow=ncol(allb1s)))
# colnames(b1Summary) <- c("site","mean","sd")
# b1Summary$site <- sites
# b1Summary$mean <- colMeans(allb1s)
# b1Summary$sd <- apply(X=allb1s,FUN=sd,MARGIN = 2)
# save(b1Summary,file="b1Summary_noTrigger.RData")
# 
# b2Summary <- data.frame(matrix(ncol=3,nrow=ncol(allb2s)))
# colnames(b2Summary) <- c("site","mean","sd")
# b2Summary$site <- sites
# b2Summary$mean <- colMeans(allb2s)
# b2Summary$sd <- apply(X=allb2s,FUN=sd,MARGIN = 2)
# save(b2Summary,file="b2Summary_noTrigger.RData")
# 
# b3Summary <- data.frame(matrix(ncol=3,nrow=ncol(allb3s)))
# colnames(b3Summary) <- c("site","mean","sd")
# b3Summary$site <- sites
# b3Summary$mean <- colMeans(allb3s)
# b3Summary$sd <- apply(X=allb3s,FUN=sd,MARGIN = 2)
# save(b3Summary,file="b3Summary_noTrigger.RData")

