############## Code executable on Jasmin #######################
# prj <- system.file("proj", package = "terra")[1]
# Sys.setenv("PROJ_LIB" = prj)

############## LIBRARIES ####################### #######################
dir_lib<-"/gws/nopw/j04/uknetzero/mesoclim/mesoclim_lib"
library(terra)
library(sf)
library(mesoclim, lib.loc=dir_lib)
library(lubridate)
# library(mesoclim)

############## 1A INPUT FILES & DIRECTORIES ####################### #######################

# basepath to badc/... oaths can be set is testing - use "" for runs on Jasmin
ceda_basepath <-""
# ceda_basepath <-"D:"

# Any plot or print outputs? set to FALSE for Jasmin runs
outputs<-FALSE
# outputs<-TRUE

# Root directory relative to these data inputs
dir_root<-"/gws/nopw/j04/uknetzero/mesoclim"
# dir_root<-"D:/"

dir_era5<-file.path(dir_root,'era5')  # ERA5 data directory

## Mac
# dir_root<-
# dir_era5<-'/Volumes/LaCie/Microclimate_2023/Data/era5'
# dir_ukcp<-'/Volumes/LaCie/Microclimate_2023/Data/ukcp18/land_rcm'
# dir_era5dly<-'/Volumes/LaCie/Microclimate_2023/Data/era5_daily'

############## 1B PARAMETERS ####################### #######################

# Timeperiod on which to carry out bias correction
pst_sdate<-as.POSIXlt('2018/01/01')
pst_edate<-as.POSIXlt('2018/12/31')

#### ADD dtmc or aoi to define area of interest !!!!!!
#dtm<-?????

############## 2 HISTORIC OBSERVATIONAL DATA ####################### #######################
### ERA5 data
# Get dtm informing era5 data
lsm<-rast(file.path(dir_era5,"lsm.tif"),subds='lsm')
geop<-rast(file.path(dir_era5,'geopotential.tif'),subds='z')
era5dtm<-mask(geop/9.80665,lsm,maskvalues=0)
plot(era5dtm)


yrs<-seq(from=year(pst_sdate),to=year(pst_edate),by=1)
# y<-2018
# Load era5 yearly data and convert to daily
for(y in yrs){
  era5file<-file.path(dir_era5,paste0('era5_surface_ukeire_',y,'.nc'))
  file.exists(era5file)
  mesoclim:::.era5todaily(filein=era5file,pathout=dir_era5dly,landsea=lsm)
  #era5toclimarray(era5file, dtm=era5dtm, aoi=dtm, dtr_cor_fac = 1.285, toArrays=TRUE)
}

### Haduk
hadukdata<-addtrees_hadukdata(pst_sdate, pst_edate,dtm, varn=c('rainfall','tasmax','tasmin'),freq='day')
# Aggregate to resolution of ukcp18 (12km)
