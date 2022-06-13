# chlorophyllCycling
This document contains the names and descriptions of the different R scrips in this repo that contain code for the modeling and analyses for the manuscript:

"A trigger may not be necessary to cause senescence in deciduous broadleaf forests" by Kathryn I Wheeler and Michael C. Dietze

Email: kwheelerecology@gmail.com

The organization of the scrips and this document as still in progress to make clearer before publications. So we apologize if things break or unclear still. 

Some of the scripts need functions from some of my R packages on github, specifically: 
library(PhenoForecast) #From github: k-wheeler/NEFI_pheno/PhenoForecast
library(PhenologyBayesModeling) #From github: k-wheeler/NEFI_pheno/PhenologyBayesModeling 
library(ecoforecastR) #From github: EcoForecast/ecoforecastR

Some of the scripts rely on reading in a csv file of site characteristics: read.csv('/projectnb/dietzelab/kiwheel/chlorophyllCycling/allPhenocamDBsitesComplete.csv',header=TRUE)
This is a csv where the columns are siteName (PhenoCam's sitename), Lat, Long, URL (URL of data archive on PhenoCam site for specific site and ROI), PFT (i.e., DB), TZ (offset from GMT, e.g., -5 for ET), startDate (first date of full year of data, e.g., "2013-01-01"), URL2 (additional PhenoCam URL if switched cameras/ROI), URL3, endDate

#Prepping Data
downloadERA5Calibration.R #Downloads ERA5 data for selected sites. You need to have the api set up to do this. 
load_ERA5.R #Files to load ERA5 met data from the saved files. 

createElmoreFitsForRescaling.R #Saves R object for each site with rescaling info as file=paste(dataDirectory,siteName,"_phenopixOutputs.RData",sep="")

createDataObjects #Saves comprehensive data objects for each site as: file=paste0(dataDirectory,siteName,"_dataFinal_includeJuly.RData"))


#Calculates PhenoCam Transition Dates
estimatePhenoCamTransitions_changePoint.R #Fits changepoint Bayesian model to PhenoCam data to estimate start of senescence transition dates and stores the model fits in the file file=paste0('varBurns/',siteName,"_",yrName,"_PhenoCam_changePointCurve_varBurn_updatedK.RData")) to be used below 
calculatePhenoCamTransitionsFromMean.R #Combines calculated PhenoCam SOS transition dates into one file: "phenocamTransitions_fromMeanFiltered.csv"

#Fit Models
createModelCalibration_climatology.R #Creates historical average null model fits 
createModelCalibrations_forecasts.R #Calibrates model to defined sites and number of included day of years

#Validation site predictions
uncertaintyAnalysisHindcasts_allSites.R #Validation site predictions for each calibration site 

#Fit Plots
plotClimatologyFits.R #Plots the historical averages model fits
plotHindcast_forecasts.R #Plots the model fits for various amounts of included data and for specified sites 
plotPhenocamTransitions.R #Plots the estimated PhenoCam transitions 

#Analyses
investigateInflectionPoints.R #Determines if inflection was met
investigateInflectionPoints_OOSsites.R #Determines if inflection was met in out of sample sites 
createCRPSpercentageMatrix.R #Creates a matrix of % site-years where our model is better than climatology over differing amounts of included data for Fig. 3c

createManuscriptFiguresFINAL.R #File to create figures for the manuscript. Figures were further edited elsewhere. 


