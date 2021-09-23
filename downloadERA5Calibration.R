library("reticulate")
library(doParallel)
n.cores <- 36

#register the cores.
registerDoParallel(cores=n.cores)
setwd("/projectnb/dietzelab/kiwheel/ERA5")

siteData <- read.csv('/projectnb/dietzelab/kiwheel/chlorophyllCycling/
                     allPhenocamDBsitesComplete.csv',header=TRUE)

cdsapi <- reticulate::import("cdsapi")
cclient <- cdsapi$Client()

end_date=as.Date("2020-12-31")

variables <- tibble::tribble(
  ~cf_name, ~units, ~api_name, ~ncdf_name,
  "air_temperature", "Kelvin", "2m_temperature", "t2m",
  "air_pressure", "Pa", "surface_pressure", NA_character_,
  NA_character_, "Kelvin", "2m_dewpoint_temperature", NA_character_,
  "precipitation_flux", "kg/m2/s", "total_precipitation", NA_character_,
  "eastward_wind", "m/s", "10m_u_component_of_wind", NA_character_,
  "northward_wind", "m/s", "10m_v_component_of_wind", NA_character_,
  "surface_downwelling_shortwave_flux_in_air", "W/m2", "surface_solar_radiation_downwards", NA_character_,
  "surface_downwelling_longwave_flux_in_air", "W/m2", "surface_thermal_radiation_downwards", NA_character_
)

var <- variables[["api_name"]][[1]]#4
foreach(i=1:nrow(siteData)) %dopar% {
  #for(i in 1:nrow(siteData)){
  #i <- 1
  siteName <- as.character(siteData$siteName[i])
  print(siteName)
  outfolder <- paste("Data/",siteName,sep="")
  # dir.create(outfolder)
  # print(paste("Created Folder:",outfolder))
  lat <- as.numeric(siteData$Lat[i])
  long <- as.numeric(siteData$Long[i])
  start_date <- as.Date(siteData$startDate[i])
  
  area <- rep(round(c(lat, long) * 4) / 4, 2)

  fname <- file.path(outfolder, paste(siteName,"_",start_date,"_",end_date,"_era5TemperatureMembers.nc", sep =""))
  #fname <- file.path(outfolder, paste(siteName,"_",start_date,"_",end_date,"_era5PrecipitationMembers.nc", sep =""))
  if(!file.exists(fname)){
    do_next <- tryCatch({
      cclient$retrieve(
        "reanalysis-era5-single-levels",
        list(
          variable = var,
          product_type = 'ensemble_members',
          date = paste(start_date, end_date, sep = "/"),
          time = "00/to/23/by/1",
          area = area,
          grid = c(0.25, 0.25),
          format = "netcdf"
        ),
        fname
      )
      FALSE
    }, error = function(e) {
      print("Failed to download variable Mean")
      TRUE
    })
  }
}


