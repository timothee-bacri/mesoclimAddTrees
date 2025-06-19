library(ecmwfr)
library(mesoclimAddTrees)
dir_out<-"D:/era5"
cds_key<-"664bf7b4-9bc2-4ce5-b479-6b731c7a4e31"
wf_set_key(key=cds_key,user='J.Mosedale@exeter.ac.uk')

for(yr in c(2024)){
  print(yr)
  download_era5(dir_out, file_out=paste0("era5_surface_ukeire_",yr,".nc"),
                          yr, area = c(61, -11.25, 48.75, 2.5),
                          variables =  c("10m_u_component_of_wind",
                                         "10m_v_component_of_wind",
                                         "2m_dewpoint_temperature",
                                         "2m_temperature",
                                         "mean_sea_level_pressure",
                                         "surface_pressure",
                                         "total_precipitation",
                                         "surface_solar_radiation_downwards",
                                         "total_sky_direct_solar_radiation_at_surface",
                                         "mean_surface_downward_long_wave_radiation_flux"),
                                         era5_user ='J.Mosedale@exeter.ac.uk',
                era5_key = "664bf7b4-9bc2-4ce5-b479-6b731c7a4e31"
                )
}

download_haduk(dir_haduk,  as.POSIXlt("2017/01/01"),as.POSIXlt("2017/12/31"),
               c("tasmax","tasmin"),freq = "day", cedausr, cedapwd, access_token
)
download_haduk(dir_haduk,  as.POSIXlt("2017/01/01"),as.POSIXlt("2017/12/31"),
               "rainfall",freq = "day", cedausr, cedapwd, access_token
)
