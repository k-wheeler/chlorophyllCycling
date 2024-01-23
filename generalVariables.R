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