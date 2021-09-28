##' Calculates cummulative Tair (within one year for one ensemble)
##'
##' @param ERA5dataFolder
##' @param calFileName
##' @param TZ_offset
##' @param variable Tair or precip
##' @import xts
##' @import ncdf4
##' @export
load_ERA5 <- function(ERA5dataFolder,calFileName,TZ_offset,variable) {
  ensembleFile <- nc_open(paste0(ERA5dataFolder,calFileName))
  if(variable=="Tair"){
    Tairs <- ncvar_get(ensembleFile)-273 #Convert from Kelvin to C
    if(!is.na(dim(Tairs)[3])){
      Tairs <- Tairs[,1,]
    }
    
    timeHours <- ensembleFile$dim$time$vals #Hours since 1900-01-01 00:00:00.0
    
    ##Convert times to actual times
    times <- as.POSIXct(timeHours*3600, origin = "1900-01-01",tz = "GMT")
    times <- times + TZ * 60 * 60 #(add number of offset hours in seconds)
    
    #Daily average
    allDates <- lubridate::date(times)
    dates <- seq(lubridate::date(times[5]),lubridate::date(times[length(times)]),"day")
    TairsDaily <- matrix(nrow=10,ncol=length(dates))
    for(d in 1:length(dates)){
      subTairs <- Tairs[,allDates==dates[d]]
      TairsDaily[,d] <- apply(subTairs,MARGIN=1,mean)
    }
    output <- TairsDaily
    nc_close(ensembleFile)
  }
  return(output)
}

