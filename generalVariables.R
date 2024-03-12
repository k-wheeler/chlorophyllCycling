headFilePath <-  '/projectnb/dietzelab/kiwheel/chlorophyllCycling/'
modelOutputFolder <- "modelOutputs/"
climatologyModelOutputsFolder <- paste0(modelOutputFolder,"climatologyModelOutputs/")
CCmodelOutputsFolder <- paste0(modelOutputFolder,'CCmodelOutputs/')
triggerModelOutputsFolder <- paste0(modelOutputFolder,'triggerModelOutputs/')
transitionEstimateOutputsFolder <- paste0(modelOutputFolder,'transitionEstimateOutputs/')

siteData <- read.csv(paste0(headFilePath,'allPhenocamDBsitesComplete.csv'),header=TRUE)
siteDataSummary <- read.csv(paste0(headFilePath,'allPhenocamDBsitesSummary.csv'),header=TRUE) #Table S1

badSiteYears <- read.csv('phenocamSitesBadYears.csv',header = FALSE) #Determined via visual assessments 
dataDirectory <- "data/"

sites <- c("harvard","umichbiological","bostoncommon","coweeta","howland2",
           "morganmonroe","missouriozarks","queens","dukehw","lacclair","bbc1","NEON.D08.DELA.DP1.00033",
           "bartlettir","oakridge1","alligatorriver","readingma",
           "bullshoals",'russellsage','sanford','downerwoods','willowcreek','laurentides','boundarywaters')
allSites <- as.character(siteData$siteName)

yearsRemoved <- c(2015,2010,2020,2015,2017,2015,2015,2010,2017,2019,2018,2017,
                  2012,2019,2019,2010,2014,2015,2017,2018,2016,2011,2012,2019)

ERA5dataFolder <- "/projectnb/dietzelab/kiwheel/ERA5/Data/"

allPhenoTranFile <- "phenocamTransitions_fromMeanFiltered.csv"
crpsMatFile <- "crpsMat_includedVsDay_Complete.RData"
reorganizedDatFile <- "reorganizedDat_includedVsDay.RData"
daysOffsetFile <- "daysOffset_includedVsDay.RData"

n.cores <- 8 #Number of cores available for parallelization
nchain=5
cutOff <- -0.02 #Cut-off For the value of second-differences that constitutes an inflection based off of PhenoCam data
NT <- 184
avgN <- 15
days <- seq(1,NT)
Nmc <- 5000
