library('rjags')
library('runjags')
source('generalVariables.R')
options(stringsAsFactors = FALSE)

output <- matrix(nrow=length(sites),ncol=7)
for(s in seq_along(sites)){
  #for(s in 1:1){
  siteName <- sites[s]
  print(siteName)
  load(paste0(dataDirectory,siteName,"_dataFinal.RData"))
  n=183
  
  fileName <- paste0(CCmodelOutputsFolder,siteName,"_",n,"_ccModel_forecast_calibration_varBurn.RData")
  load(fileName) #Load Model Output 
  
  out.mat <- data.frame(as.matrix(out.burn$params))
  
  b0 <- out.mat$b0
  b4 <- -1*out.mat$b4 #Model in paper has a negative where this model was calibrated with positive
  b3 <- out.mat$b3
  output[s,] <- c(siteName,round(mean(b0),digits=5),
                  round(sd(b0),digits=5),
                  round(mean(b3),digits=5),
                  round(sd(b3),digits=5),
                  round(mean(b4),digits=5),
                  round(sd(b4),digits=5))

}

colnames(output) <- c("Site Name","b0_mean","b0_sd","b1_mean","b1_sd",
                      "b2_mean","b2_sd")   

output <- output[order(as.numeric(output[,4])),]
write.csv(output,file="modelParameterCalibrations.csv",quote = FALSE,row.names = FALSE)
