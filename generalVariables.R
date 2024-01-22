headFilePath <-  '/projectnb/dietzelab/kiwheel/chlorophyllCycling/'
modelOutputFolder <- "modelOutputs/"
climatologyModelOutputsFolder <- paste0(modelOutputFolder,"climatologyModelOutputs/")
CCmodelOutputsFolder <- paste0(modelOutputFolder,'CCmodelOutputs/')
transitionEstimateOutputsFolder <- paste0(modelOutputFolder,'transitionEstimateOutputs/')

siteData <- read.csv(paste0(headFilePath,'allPhenocamDBsitesComplete.csv'),header=TRUE)

badSiteYears <- read.csv('phenocamSitesBadYears.csv',header = FALSE) #Determined via visual assessments 
dataDirectory <- "data/"

sites <- c("harvard","umichbiological","bostoncommon","coweeta","howland2",
           "morganmonroe","missouriozarks","queens","dukehw","lacclair","bbc1","NEON.D08.DELA.DP1.00033",
           "arbutuslake","bartlettir","proctor","oakridge1","hubbardbrook","canadaOA","alligatorriver","readingma",
           "bullshoals","thompsonfarm2N","ashburnham","shalehillsczo")

yearsRemoved <- c(2015,2010,2020,2015,2017,2015,2015,2010,2017,2019,2018,2017,
                  2012,2019,2019,2010,2014,2015,2017,2018,2016,2011,2012,2019)

ERAdataFolder <- "/projectnb/dietzelab/kiwheel/ERA5/"

allPhenoTranFile <- "phenocamTransitions_fromMeanFiltered.csv"
crpsMatFile <- "crpsMat_includedVsDay_Complete.RData"
reorganizedDatFile <- "reorganizedDat_includedVsDay.RData"
daysOffsetFile <- "daysOffset_includedVsDay.RData"

n.cores <- 8 #Number of cores available for parallelization
nchain=5

#General File Names:
#ERA5 data name: paste0(ERAdataFolder,siteName,"_",start_date,"_",end_date,"_era5TemperatureMembers.nc")
#Rescaled parameters: file=paste0(dataDirectory,siteName,"_phenopixOutputs.RData"))
#Data Objects: file=(paste0(dataDirectory,siteName,"_dataFinal.RData"))
#Changepoint model outputs to estimate date of inflection: paste0(transitionEstimateOutputsFolder,siteName,"_",yrName,"_PhenoCam_changePointCurve_varBurn.RData")
#Climatology model outputs: paste0(climatologyModelOutputsFolder,siteName,"_climatology_forecast_calibration_varBurn.RData")
#CC model outputs: paste0(CCmodelOutputsFolder,siteName,"_",n,"_ccModel_forecast_calibration_varBurn.RData")
#CRPS values of out of sample predictions:paste0("outOfSampleSites_crps_",calSite,"_183.csv")
#Inflection for calibration sites: paste0(siteName,"_inflectionPointData_15.RData"))
#Inflection for out of sample predictions: paste0(calSite,"_inflectionPointData_OOSsites_mean15_",n,".RData")
