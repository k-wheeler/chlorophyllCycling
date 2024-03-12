calculateStart <- function(ys,avgNum=15){
  avgDiffs <- rep(0,avgNum)
  for(t in (avgNum + 1):(length(ys)-avgNum)){
    prevAvg <- mean(diff(ys[1:(t-1)]),na.rm = TRUE)
    newAvg <- mean(diff(ys[t:(t+avgNum)]),na.rm = TRUE)
    avgDiffs <- c(avgDiffs,(newAvg-prevAvg))
  }
  avgDiffs <- c(avgDiffs,rep(0,avgNum))
  return(avgDiffs)
}