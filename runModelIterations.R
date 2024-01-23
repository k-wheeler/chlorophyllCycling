##' @name mat2mcmc.list
##' @title mat2mcmc.list
##' @export
##' @author Mike Dietze
##' @description convert a matrix to a CODA mcmc.list
##' @param w  matrix
mat2mcmc.list <- function(w) {
  temp <- list()
  chain.col <- which(colnames(w) == "CHAIN")
  for (i in unique(w[, "CHAIN"])) {
    temp[[i]] <- coda:::as.mcmc(w[w[, "CHAIN"] == i, -chain.col])
  }
  return(as.mcmc.list(temp))
}

##' Adds additional iterations until convergence and large enough effective sample size are met
##'
##' @param j.model The JAGS model
##' @param variableNames Variable names for coda.samples
##' @param maxIter The maximum number of iterations to run before declaring that the model will not converge (optional)
##' @param baseNum The initial number of iterations to run before checking for convergence
##' @param iterSize The number of iterations to iteratively increase before checking for convergence again and large enough sample size (optional)
##' @param ID A identification string to paste with the message of non-convergence (optional)
##' @param effSize The desired necessary effective sample size (default is 5000 samples)
##' @param partialFile The desired name to save partial output (non-converged) to
##' @param startFileName The file name of a previous jags output that you want to add to (used to combine new runs with previous samples) (Note that this variable must be "partialOutput")
##' @export
##' @import rjags
##' @import runjags
##' @import coda
runForecastIter <- function(j.model,variableNames,maxIter=10**9,baseNum=5000,iterSize =5000,ID="",effSize=5000,partialFile=FALSE,startFileName=FALSE){
  if(typeof(startFileName)==typeof(FALSE)){
    jags.out   <- coda.samples (model = j.model,
                                variable.names = variableNames,
                                n.iter = baseNum)
    
    numb <- baseNum #The current number of iterations at this step
  }else{
    load(startFileName)
    #loaded as partialOutput
    baseNum <- nrow(partialOutput)
    partialOutputCombined <- cbind(as.matrix(partialOutput$params),as.matrix(partialOutput$predict,chains=TRUE))
    jags.out <- mat2mcmc.list(partialOutputCombined)
  }
  
  ###Check for Model Convergence and keep adding to MCMC chains if it hasn't converged and/or effective sample size is not large enough
  
  continue <- TRUE #Flag that signals that the coda.samples should be rerun to produce more iterations because the model hasn't converged yet/doesnt have a large enough sample size
  GBR.bad <- TRUE #Flag that signals that it hasn't converged yet
  burnin <- 0
  while(continue & numb<maxIter){
    print(numb) #Just done to keep track of the number of iterations that have been performed
    new.out   <- coda.samples(model = j.model,
                               variable.names = variableNames,
                               n.iter = iterSize)
    numb <- numb + iterSize
    jags.out <- combine.mcmc(mcmc.objects=list(jags.out,new.out),collapse.chains = FALSE)
    continue <- FALSE
    if(GBR.bad){
      out = list(params=NULL,predict=NULL) #Split output into parameters and state variables
      mfit = as.matrix(jags.out,chains=TRUE)
      pred.cols = grep("x[",colnames(mfit),fixed=TRUE)
      chain.col = which(colnames(mfit)=="CHAIN")
      out$params = mat2mcmc.list(mfit[,-pred.cols])
      GBR.vals <- gelman.diag(out$params)
      print(GBR.vals$psrf)
      for(i in 1:nrow(GBR.vals$psrf)){ #Checking to see if any of the parameters haven't converged yet
        for(j in 1:ncol(GBR.vals$psrf)){
          if(!is.nan(GBR.vals$psrf[i,j])){
            if(GBR.vals$psrf[i,j]>1.04){
              continue <-  TRUE
            }
          }
        }
      }
      if(typeof(partialFile)!=typeof(FALSE)){
        print("entered if")
        partialOutput <- list(params=out$params)
        partialOutput$predict <- mat2mcmc.list(mfit[,c(chain.col,pred.cols)])
        
        out.mat <- as.matrix(partialOutput$params)
        thinAmount <- round(nrow(out.mat)/5000,digits=0)
        partialOutput2 <- list()
        partialOutput2$params <- window(partialOutput$params,thin=thinAmount)
        partialOutput2$predict <- window(partialOutput$predict,thin=thinAmount)
        partialOutput <- partialOutput2
        save(file=partialFile,partialOutput)
        print("saved output")
      }
      
    }
    if(!continue){
      if(burnin==0){ #If the while loop has to be rerun because the effective size is too small, you don't need to calculate burnin again
        GBR <- gelman.plot(out$params)
        burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[,,2]>1.05,1,any)),1)+1]
        if(length(burnin) == 0) burnin = 1
      }
      var.burn <- window(jags.out,start=burnin)
      out.burn = list(params=NULL,predict=NULL)
      mfit = as.matrix(var.burn,chains=TRUE)
      pred.cols = grep("x[",colnames(mfit),fixed=TRUE)
      chain.col = which(colnames(mfit)=="CHAIN")
      out.burn$params = mat2mcmc.list(mfit[,-pred.cols])
      effsize <- effectiveSize(out.burn$params)
      print(effsize)
      for(i in 1:length(effsize)){
        if(effsize[i]<effSize){
          continue = TRUE
        }
      }
    }
  }
  if(!continue){
    out.burn$predict = mat2mcmc.list(mfit[,c(chain.col,pred.cols)])
    return(out.burn)
  }
  else{
    print(paste(ID, "model did not converge, returning FALSE"))
    return(FALSE)
  }
}

##' Continues to run MCMC until convergence and large enough effective sample size
##'
##' @param j.model The jags model
##' @param variableNames A vector of the variable names that need to converge
##' @param maxIter The maximum number of iterations to run
##' @param baseNum The number of initial iterations to run
##' @param iterSize The number of iterations per each run (between when GBR values are checked)
##' @param maxGBR The maximum allowable GBR value after the baseNum of iterations are run
##' @param ID An optional identifier to print if it doesn't converge
##' @param sampleCutoff The minimum desired effective sample size for testing (default is 5000)
##' @import rjags
##' @import runjags
##' @import coda
##' @export
runMCMC_Model <- function(j.model,variableNames,maxIter=1000000000,baseNum=80000,iterSize=40000,maxGBR=10,ID="",sampleCutoff=5000){
  var.out   <- coda.samples(model = j.model,
                             variable.names = variableNames,
                             n.iter = baseNum)
  numb <- baseNum
  continue <- TRUE
  GBR.bad <- TRUE
  burnin <- 0
  while(continue & numb<maxIter){
    print(numb)
    new.out   <- coda.samples(model = j.model,
                              variable.names = variableNames,
                              n.iter = iterSize)
    var.out <- combine.mcmc(mcmc.objects=list(var.out,new.out),collapse.chains = FALSE)
    continue <- FALSE
    if(GBR.bad){
      GBR.vals <- gelman.diag(var.out)
      GBR.bad <- FALSE
      for(i in 1:nrow(GBR.vals$psrf)){
        for(j in 1:ncol(GBR.vals$psrf)){
          if(GBR.vals$psrf[i,j]>maxGBR){
            print(GBR.vals)
            print(c("GBR values too high:",ID))
            return(FALSE)
          }
          if(GBR.vals$psrf[i,j]>1.04){
            continue <-  TRUE
            GBR.bad <- TRUE
          }
        }
      }
      print(GBR.vals)
    }
    if(!continue){
      if(burnin==0){
        GBR <- gelman.plot(var.out)
        burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[,,2]>1.05,1,any)),1)+1]
        if(length(burnin) == 0) burnin = 1
      }
      var.burn <- window(var.out,start=burnin)
      #var.burn <- var.out
      effsize <- effectiveSize(var.burn)
      for(i in 1:length(effsize)){
        if(effsize[i]<sampleCutoff){
          continue = TRUE
        }
      }
      print(effsize)
    }
    numb <- numb+iterSize
  }
  if(continue==TRUE){
    print("Model Did not Converge")
    var.burn <- FALSE
  }
  return(var.burn)
}