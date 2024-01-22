library("reticulate")
library(doParallel)
source('generalVariables.R')
n.cores <- 24

#register the cores.
registerDoParallel(cores=n.cores)

cdsapi <- reticulate::import("cdsapi") #Need to install separately so you can interact with api
cclient <- cdsapi$Client()

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

var <- variables[["api_name"]][[1]] #Only Download Air temperature
foreach(i=1:nrow(siteData)) %dopar% {
  siteName <- as.character(siteData$siteName[i])
  print(siteName)

  lat <- as.numeric(siteData$Lat[i])
  long <- as.numeric(siteData$Long[i])
  start_date <- as.Date(siteData$startDate[i])
  end_date <- as.Date(siteData$endDate[i])
  
  area <- rep(round(c(lat, long) * 4) / 4, 2)
  
  fname <- file.path(ERAdataFolder, paste0(siteName,"_",start_date,"_",end_date,"_era5TemperatureMembers.nc"))

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
