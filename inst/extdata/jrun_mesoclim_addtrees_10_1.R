#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

############## Check User PARAMETERS passed via Rscript line command ######
# 1 parcel_file location
# 2 output directory
# 3 model memeber number eg "01"
# 4 start year
# 5 end year

if (length(args)==0) {
  stop("At least one argument must be supplied (parcel file)", call.=FALSE)
} else if (length(args)>4) {
  stop("Up to four arguments are required (parcel file, output dir, starttime, endtime)", call.=FALSE) }


if(!file.exists(args[1])) stop("Input file provided does NOT exist!!!")

maxyear<-2080
minyear<-1981
out_default<-"/gws/nopw/j04/uknetzero/mesoclim/mesoclim_outputs"

if (length(args)<5) {
  warning("Using default end times!!")
  args[5] = 2080
}
if (length(args)<4) {
  warning("Using default start and end times!!")
  args[4] = 2024
}
if (length(args)<3){
  warning("Using default UKCP modelrun (01)!!")
  args[3] = "01"

}
if (length(args)<2) {
  warning("Using default output directory!!")
  args[2] = out_default
}

if (class(args[4])[1]!="numeric") stop("Invalid 4th parameter!!")
if (class(args[5])[1]!="numeric") stop("Invalid 5th parameter!!")

if(!dir.exists(args[2])){
  warning(paste("Output directory provided does NOT exist, using default dir",out_default))
  args[2]<-out_default
}
if(args[5]>maxyear){
  warning(paste("End year too high - changing to",maxyear))
  args[5]<-maxyr
}
if(args[4]<minyear){
  warning(paste("Start year too low - changing to",minyear))
  args[4]<-minyr
}

############## Assign PARAMETERS and define others ####################### #######################
parcels_file<-args[1]
dir_out<-args[2]
modelrun<-args[3]
ftr_sdate<-as.POSIXlt(paste0(args[4],'2021/01/01'))
ftr_edate<-as.POSIXlt(paste0(args[5],'/12/31'))


# These are fixed for ADDTREES analyses - shouldn't need to change
collection<-'land-rcm'
domain<-'uk'


############## LIBRARIES ####################### #######################
dir_lib<-"/gws/nopw/j04/uknetzero/mesoclim/mesoclim_lib"
library(terra)
library(sf)
library(mesoclim, lib.loc=dir_lib)
library(lubridate)
library(mesoclimAddTrees, lib.loc=dir_lib)

terraOptions(tempdir = "/gws/nopw/j04/uknetzero/mesoclim/terra_storage")

############## 1A INPUT FILES & DIRECTORIES ####################### #######################

# basepath to badc/... oaths can be set is testing - use "" for runs on Jasmin
ceda_basepath <-""

# Any plot or print outputs? set to FALSE for Jasmin runs
outputs<-FALSE

# Root directory relative to these data inputs
dir_root<-"/gws/nopw/j04/uknetzero/mesoclim"

# Filepath to vector file of parcels output by ellicitor app.
parcels_file<-file.path(dir_root,'mesoclim_inputs','parcels','land_parcels.shp') # elicitor app output file
if(!file.exists(parcels_file)) stop("Cannot find parcels input file!!")

# Filepath to coastline boundary polygon (not necessary if dtm already masked)
coast_file<-file.path(dir_root,'mesoclim_inputs','boundaries','CTRY_DEC_2023_UK_BGC.shp') # MHW line generalised to 20m
if(!file.exists(coast_file)) stop("Cannot find coastal mask file for the UK!!")

# Filepath to fine resolution DTM of UK (OS Terrain50 -  which has advantages in having a .vrt raster that can be queried to enhance cropping and extraction of relevant area to match aoi etc. Subdirectories within hold actual tiled data.)
#ukdtm_file<-file.path(dir_root,'mesoclim_inputs','dtm',"GBdem50.vrt") # 50m dtm virtual raster
ukdtm_file<-file.path(dir_root,'mesoclim_inputs','dtm',"sw_dtm.tif") # 50m dtm virtual raster
if(!file.exists(ukdtm_file)) stop("Cannot find DTM file for the UK!!")

# Directory for outputs - to which individual parcel .csv timeseries files are written.
dir_out<-file.path(dir_root,'mesoclim_outputs')  # output dir
if(!dir.exists(dir_out)) stop("Cannot find output directory!!")

# Remove any existing parcel files in dir_out
#pfiles<-list.files(dir_out,full.names=TRUE, pattern="parcel")
#file.remove(pfiles)


############## 2 PREPARE INPUTS ####################### #######################
# TO DO: Prepare climate data for whole decade - more effecicient than opening file for every year

### Area of interest and elevation data - AOI for downscaling defined by parcel data

# Load UK fine resolution dtm & boundary(coast)
dtmuk<-terra::rast(ukdtm_file) # 50m resolution dtm of all UK
coast_v<-terra::project(terra::vect(coast_file),dtmuk)

# Load parcels file and project to crs of output dtm (OS coords)
parcels_v<-terra::project(terra::vect(parcels_file),dtmuk)

# Generate local area and dtm and wider extents
aoi<-terra::vect(terra::ext(parcels_v))
terra::crs(aoi)<-terra::crs(dtmuk)

# Load ukcp coarse resolution dtm for aoi
dtmc<-get_ukcp_dtm(aoi, basepath=ceda_basepath)

# Create fine resolution dtm of downscaling area  - ensure they fall within extent of loaded dtm & mask to coast_v (sets sea to NA)
dtmf<-terra::mask(terra::crop(terra::crop(dtmuk,aoi),dtmuk),coast_v)

# Generate medium area and resoilution dtm (for coatal/wind effects)
dtmm<-get_dtmm(aoi,dtmc,dtmuk,basepath=ceda_basepath)

# Plot dtmf and overlay parcels
if(outputs){
  plot(dtmc,main='DTMs')
  plot(dtmm,add=TRUE)
  plot(dtmf,add=TRUE)
  plot(parcels_v,add=TRUE)
}

### Prepare climate and seas surface data
t0<-now()

# Process climate data from UKCP18 regional files on ceda archive
climdata<-addtrees_climdata(aoi,ftr_sdate,ftr_edate,collection='land-rcm',domain='uk',member='01',basepath=ceda_basepath)

# Process sea surface temo data from ceda archive ??CHANGE OUTPUT TO PROJECTION OF DTMC/AOI? COMBINE WITH ABOVE?
sstdata<-addtrees_sstdata(ftr_sdate,ftr_edate,aoi=climdata$dtm,member='01',basepath=ceda_basepath)

dataprep_time<-now()-t0
print(paste("Time for preparing data =", format(dataprep_time)))


if(outputs){
  plot(project(sstdata[[1]],crs(dtmc)))
  plot(dtmm,add=T)
  plot(aoi,add=TRUE)
}

# Check data - plot summary figs
if(outputs) climdata<-checkinputs(climdata, tstep = "day")


############## 3 SPATIAL DOWNSCALE ####################### #######################
years<-unique(c(year(ftr_sdate):year(ftr_edate)))
for (yr in years){
# To DO: Add yearly loop - downscale->parcelcalcs
  sdatetime<-as.POSIXlt(paste0(yr,'/01/01'))
  edatetime<-as.POSIXlt(paste0(yr,'/12/31'))
  mesoclimate<-spatialdownscale(subset_climdata(climdata,sdatetime,edatetime), subset_climdata(sstdata,sdatetime,edatetime),
                                dtmf, dtmm, basins = NA, cad = TRUE,
                                coastal = TRUE, thgto =2, whgto=2,
                                rhmin = 20, pksealevel = TRUE, patchsim = TRUE,
                                terrainshade = FALSE, precipmethod = "Elev", fast = TRUE, noraincut = 0.01)
  downscale_time<-now()-t0
  print(paste("Time for downscaling =", format(downscale_time)))

  if(outputs){
    climvars<-c('tmin','tmax','relhum','pres','swrad','lwrad','windspeed','winddir','prec')
    for(var in climvars){
      print(var)
      r<-mesoclimate[[var]]
      names(r)<-rep(var,nlyr(r))
      plot_q_layers(r,vtext=var)
    }
  }

# write_climdata(mesoclimate,file.path(dir_out,'mesoclimate_1yr_test.Rds'))


##############  4 Calculate and write parcel outputs  ####################### #######################
#Calculate weighted means of climate variables for each parcel and write as .csv files.

# Calculate parcel values
#t0<-now()
  parcel_list<-create_parcel_list(mesoclimate,parcels_v,id='gid')
  write_parcels(parcel_list, dir_out, overwrite='append')
  parcel_time<-now()-t0
  print(paste("Time for parcel calculation and writing =", format(parcel_time)))
}
total_time<-now()-t0
print(paste("Total time =", format(parcel_time)))

#  var_sf<-get_parcel_var(mesoclimate,'tmax', parcels_v,id='gid', stat='mean' )
#  map_parcel_var(var_sf, plotvar='tmax', idvar='gid')
