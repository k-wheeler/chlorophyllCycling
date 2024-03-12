# chlorophyllCycling
This document contains the names and descriptions of the different R scrips in this repo that contain code for the modeling and analyses for the manuscript:

"A trigger may not be necessary to cause senescence in deciduous broadleaf forests" by Kathryn I Wheeler and Michael C. Dietze

Email: kwheelerecology@gmail.com

The scripts require a site data csv file ('allPhenocamDBsitesComplete.csv') where the columns are siteName (PhenoCam's sitename), Lat, Long, URL (URL of data archive on PhenoCam site for specific site and ROI), PFT (i.e., DB), TZ (offset from GMT, e.g., -5 for ET), startDate (first date of full year of data, e.g., "2013-01-01"), URL2 (additional PhenoCam URL if switched cameras/ROI), URL3, endDate

The Phenocam site-years that were visually assessed as poor quality due to large gaps of missing gaps are in the file 'phenocamSitesBadYears.csv'

In order to run code on your computer, file paths should be changed in the generalVariables.R file. This file also has other general variables and is sourced in all other R scripts. 

Note: In this code, the parameter values of b0, b3, and b4 correspond to b0, b1, and -1*b2 in the manuscript


#Prepping Data
* downloadERA5Calibration.R #Downloads ERA5 data for selected sites. You need to have the api set up to do this. 

* downloadPhenocam.R #Downloaded phenocam gcc data based off of URLS in site data file

* load_ERA5.R #Files to load ERA5 met data from the saved files. 

* createElmoreFitsForRescaling.R #Saves R object for each site with rescaling info as file=paste0(dataDirectory,siteName,"_phenopixOutputs.RData"))

* createDataObjects #Saves comprehensive data objects for each site as: file=(paste0(dataDirectory,siteName,"_dataFinal.RData"))


#Calculates PhenoCam Transition Dates
* estimatePhenoCamTransitions_changePoint.R #Fits changepoint Bayesian model to PhenoCam data to estimate start of senescence transition dates and stores the model fits in the file = paste0(transitionEstimateOutputsFolder,siteName,"_",yrName,"_PhenoCam_changePointCurve_varBurn.RData") to be used below 

* calculatePhenoCamTransitionsFromMean.R #Combines calculated PhenoCam SOS transition dates into one file: "phenocamTransitions_fromMeanFiltered.csv"


#Fit Models
* runModelIterations.R #Includes two functions to iteratively run below models until convergence is reached with large enough effective sample sizes

* createModelCalibration_climatology.R #Creates historical average null model fits that are saved in file=paste0(climatologyModelOutputsFolder,siteName,"_climatology_forecast_calibration_varBurn.RData")

* createModelCalibrations_CCmodel.R #Calibrates chlorophyll cycling model to defined sites and number of included day of years and saves it each in file=paste0(CCmodelOutputsFolder,siteName,"_",n,"_ccModel_forecast_calibration_varBurn.RData")

* createModelCalibrations_triggerModel.R #STILL IN DEVELOPMENT


#Validation site predictions
* uncertaintyAnalysisHindcasts_allSites.R #Validation site predictions for each calibration site 

#Fit Plots
* plotClimatologyFits.R #Plots the historical averages model fits

* plotHindcast_forecasts.R #Plots the model fits for various amounts of included data and for specified sites 

* plotPhenocamTransitions.R #Plots the estimated PhenoCam transitions 

#Analyses
* investigateInflectionPoints.R #Determines if inflection was met and saves each in file=paste0(siteName,"_inflectionPointData_15.RData"))

* investigateInflectionPoints_OOSsites.R #Determines if inflection was met in out of sample sites and saves each in file=paste0(calSite,"_inflectionPointData_OOSsites_mean15_",n,".RData")

* createCRPSpercentageMatrix.R #Creates a matrix of % site-years where our model is better than climatology over differing amounts of included data for Fig. 3c and creates files "crpsMat_includedVsDay_Complete.RData", "reorganizedDat_includedVsDay.RData", and "daysOffset_includedVsDay.RData"

* testingImportanceOfTiming.R #Code to run a random forest model to investigate if the environmental conditions after the transition are a better predictor of the transition than those before


#Create Figures
* createManuscriptFiguresFINAL.R #File to create figures for the manuscript. Figures were further edited elsewhere. 

* createCalibrationParameterTable.R #Creates supplementary table of parameter values for each calibration site


#Extra
* ciEnvelope.R #Function to plot credible interval polygon on figures


The file naming structure is below: 

* ERA5 data name: paste0(ERAdataFolder,siteName,"_",start_date,"_",end_date,"_era5TemperatureMembers.nc")

* Rescaled parameters: file=paste0(dataDirectory,siteName,"_phenopixOutputs.RData"))

* Data Objects: file=(paste0(dataDirectory,siteName,"_dataFinal.RData"))

* Changepoint model outputs to estimate date of inflection: paste0(transitionEstimateOutputsFolder,siteName,"_",yrName,"_PhenoCam_changePointCurve_varBurn.RData")

* Climatology model outputs: paste0(climatologyModelOutputsFolder,siteName,"_climatology_forecast_calibration_varBurn.RData")

* CC model outputs: paste0(CCmodelOutputsFolder,siteName,"_",n,"_ccModel_forecast_calibration_varBurn.RData")

* CRPS values of out of sample predictions:paste0("outOfSampleSites_crps_",calSite,"_183.csv")

* Inflection for calibration sites (For each site, year, and number of days per year included combination model run, if the model converged well and what was the lowest second difference in the predictions to indicate curvature): paste0(siteName,"_inflectionPointData_15.RData"))

* Inflection for out of sample predictions: paste0(calSite,"_inflectionPointData_OOSsites_mean15_",n,".RData")
