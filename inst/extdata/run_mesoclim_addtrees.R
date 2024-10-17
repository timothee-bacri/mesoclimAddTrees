############## Code executable on LOCAL machine with downloaded inputs #######################

############## LIBRARIES ####################### #######################
# install_gitgub("ilyamaclean/mesoclim")
# install_gitgub("jrmosedale/mesoclimAddTrees")
library(terra)
library(sf)
library(mesoclimAddTrees)
library(lubridate)
library(mesoclim)

############## 1A INPUT FILES & DIRECTORIES ####################### #######################

# Any plot or print outputs? set to FALSE for Jasmin runs/slightly faster runs
#outputs<-FALSE
outputs<-TRUE

# Root directory - other dirs relative to this one or just use absolute paths below
dir_root<-"D:"

# Filepath to vector file of parcels (produced by ellicitor app - defines AOI and Outputs)
parcels_file<-file.path(dir_root,'mesoclim_inputs','parcels','land_parcels.shp') # elicitor app output file
# parcels_file<-file.path(dir_root,'mesoclim_inputs','parcels','PDNPA_Boundary.shp')
file.exists(parcels_file)

# Filepath to UK coastline boundary polygon (not necessary if dtm already masked)
coast_file<-file.path(dir_root,'mesoclim_inputs','boundaries','CTRY_DEC_2023_UK_BGC.shp') # MHW line generalised to 20m
file.exists(coast_file)

# Filepath to fine resolution DTM of UK (OS Terrain50)
ukdtm_file<-file.path(dir_root,'mesoclim_inputs','dtm',"uk_dtm.tif") # 50m dtm  raster
file.exists(ukdtm_file)

# Filepath to UKCP18 orography (coarse resolution DTM) file as downloaded from ceda
ukcpdtm_file<-file.path(dir_root,"mesoclim_inputs","ukcp18rcm","orog_land-rcm_uk_12km_osgb.nc")
file.exists(ukcpdtm_file)

# Directory for UKCP18 RCM input files
dir_ukcp<-file.path(dir_root,"mesoclim_inputs","ukcp18rcm")
dir.exists(dir_ukcp)

# Directory for UKCP18 Sea Surface temp input files
dir_sst<-file.path(dir_root,"mesoclim_inputs","ukcp18sst")
dir.exists(dir_sst)

# Directory for OUTPUTS - to which individual parcel .csv timeseries files are written.
dir_out<-file.path(dir_root,'mesoclim_outputs')  # output dir
dir.exists(dir_out)

############## 1B OTHER PARAMETERS ####################### #######################

# Start time for future climate timeseries.
startdate<-as.POSIXlt('2022/01/01')

# End time for future climate timeseries.
enddate<-as.POSIXlt('2022/12/31') # If using shared data folder use max value of as.POSIXlt('2039/12/31')

# Model run of UKCP18rcm to be downscaled.
member<-c('01')

# These are fixed for ADDTREES analyses - shouldn't need to change
collection<-'land-rcm'
domain<-'uk'


############## 2 PREPARE INPUTS ####################### #######################

###### Area of interest and elevation data - AOI for downscaling defined by parcel data

# Load UK fine resolution dtm & boundary(coast)
dtmuk<-terra::rast(ukdtm_file) # 50m resolution dtm of all UK
coast_v<-terra::project(terra::vect(coast_file),dtmuk)

# Load parcels file and project to crs of output dtm (OS coords)
parcels_v<-terra::project(terra::vect(parcels_file),dtmuk)
#  parcels_v<-crop(parcels_v,c(392883.1,432454.3,385000,411258.4)) # To crop to about 1/3 of Peak district area

# Generate local area and dtm and wider extents
aoi<-terra::vect(terra::ext(parcels_v))
terra::crs(aoi)<-terra::crs(parcels_v)

# Load ukcp coarse resolution dtm for aoi
dtmc<-get_ukcp_dtm(aoi, ukcpdtm_file)

# Create fine resolution dtm of downscaling area  - ensure they fall within extent of loaded dtm & mask to coast_v (sets sea to NA)
dtmf<-terra::mask(terra::crop(terra::crop(dtmuk,aoi),dtmuk),coast_v)

# Generate medium area and resoilution dtm (for coatal/wind effects)
dtmm<-get_dtmm(dtmf,dtmc,dtmuk)

# Plot dtmf and overlay parcels
if(outputs){
  plot(dtmc,main='DTMs')
  plot(dtmm,add=TRUE)
  plot(dtmf,add=TRUE)
  plot(parcels_v,add=TRUE)
}

# Size of area and number of parcels that it will downscale for
print(paste("Downscaling for an area of",round(expanse(aoi,"km"),0),"km^2, generating output for",nrow(parcels_v),"parcels"))


###### Calculate topographical properties - worth it if looping over several downscaling calls (eg multiple years)

# Windshelter coef
t0<-now()
wca<-calculate_windcoeffs(dtmc,dtmm,dtmf,zo=2)

# Cold air drainage basins - as above and ONLY if using coastal correction - can take several minutes for large areas
basins<-basindelin(dtmf, boundary = 2)
print(now()-t0)


####### Prepare climate and seas surface data - prepare for whole time period if up to ~ 10 years
t0<-now()

# Process climate data from UKCP18 regional files on ceda archive - PROVIDE dtmc AND aoi as ONE
climdata<-ukcp18toclimarray(dir_ukcp,dtmc,startdate,enddate,collection,domain,member,temp_hgt=2, wind_hgt=2,)

# Process sea surface temo data from ceda archive ??CHANGE OUTPUT TO PROJECTION OF DTMC/AOI? COMBINE WITH ABOVE?
sstdata<-create_ukcpsst_data(dir_sst,startdate,enddate,dtmc,member)

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
# If modelling many years will require looping for each year - appends new data to parcel file outputs with each loop
years<-unique(c(year(startdate):year(enddate)))
#yr<-years
for (yr in years){
  sdatetime<-as.POSIXlt(paste0(yr,'/01/01'))
  edatetime<-as.POSIXlt(paste0(yr,'/12/31'))

  # If cad = TRUE can take long time to calculate cold air drainage across large areas !!!
  t0<-now()
  mesoclimate<-spatialdownscale(subset_climdata(climdata,sdatetime,edatetime), subset_climdata(sstdata,sdatetime,edatetime),
                                dtmf, dtmm, basins = basins, wca=wca, cad = TRUE,
                                coastal = TRUE, thgto =2, whgto=2,
                                rhmin = 20, pksealevel = TRUE, patchsim = TRUE,
                                terrainshade = FALSE, precipmethod = "Elev", fast = TRUE, noraincut = 0.01)

  downscale_time<-now()-t0
  print(paste("Time for downscaling single year =", format(downscale_time)))

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

  ###### Calculate and write parcel outputs with each yearly iteration

  # Calculate parcel values
  t0<-now()
  parcel_list<-create_parcel_list(mesoclimate,parcels_v,id='gid')
  write_parcels(parcel_list, dir_out, overwrite='append')
  parcel_time<-now()-t0
  print(paste("Time for parcel calculation and writing =", format(parcel_time)))
  total_time<-now()-t0
  print(paste("Total parcel writing time =", format(parcel_time)))

} # end year loop


############## Get a parcel variable and plot a map of it   ####################### #######################
var_sf<-get_parcel_var(mesoclimate,'tmax', parcels_v,id='gid', stat='mean' )
map_parcel_var(var_sf, plotvar='tmax', idvar='gid')


