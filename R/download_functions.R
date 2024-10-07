# ------------------------------------ ------------------------------------ ------------------
# ------------------------------------ UKCP & HADUK DOWNLOADS ------------------------------------
# ------------------------------------ ------------------------------------ ------------------

#' @title Downloads 1km gridded historic climata data for the UK
#' @description The function `download_hadukdaily` downloads 1km HadUK-Grid Gridded
#' Climate Observation data on a 1km grid over the UK from CEDA.
#' @param dir_out directory to which to save data
#' @param startdate - POSIXlt value indicating start date for data required
#' @param enddate - POSIXlt value indicating end date for data required
#' @param varn vector of variables required, one or more of: 'rainfall','tasmax','tasmin','hurs','psl','pv','sun','sfsWind'.
#'  NOTE: Only `rainfall` (precipitation in mm), `tasmax`(maximum daily temperature, deg C) and
#' `tasmin` (minimum daily temperature, deg C) are available at daily time steps (all other monthly only).
#' @param freq - string of frequency required - either 'day' or 'mon' (monthly)
#' @param cedausr ceda username string
#' @param cedapwd ceda password sting
#' @param access_token ceda access token as provided by \href{https://services-beta.ceda.ac.uk/account/token/?_ga=2.221117848.1874942177.1727205775-405001952.1710779653}{Ceda access token generator}
#' @return names of downloaded files
#' @details Uses ceda dap service.
#' @import curl
#' @export
#' @keywords download
#' @examples
#' \dontrun{
#' dir_out <- tempdir()
#' cedausr<-"your_user_name"
#' cedapwd <- "your_password"
#' access_token<-"your_token"
#' download_haduk(dir_out,as.POSIXlt('2018-01-01'),as.POSIXlt('2018-12-31'),varn=c('tasmax','tasmin'),freq='day', cedausr,cedapwd, access_token)
#' }
download_haduk<-function(dir_out, startdate, enddate,
                         varn=c('rainfall','tasmax','tasmin'),
                         freq=c('day','mon'),
                         cedausr,cedapwd,access_token
) {
  # Checks
  varn<-match.arg(varn, several.ok=TRUE)
  freq<-match.arg(freq)
  if (any(!varn %in% c("rainfall","tasmax","tasmin")) & freq=='day') stop("Daily data are not available for chosen variables!!" )
  if(class(startdate)[1]!="POSIXlt" | class(enddate)[1]!="POSIXlt") stop("Date parameters NOT POSIXlt class!!")

  # Derive months of data required
  dateseq<-seq(startdate,enddate,'month')
  yrs<-year(dateseq)
  mnths<-month(dateseq)

  # Check time period requested
  tme<-as.POSIXlt(Sys.time())
  current_yr<-tme$year+1900
  if (any(yrs > current_yr)) stop("Function downloads observed data. Cannot be used for future years")
  if (any(yrs == current_yr)) warning("Data may not yet be available for current year")
  if ("rainfall" %in% varn) {
    if (any(yrs < 1891)) stop("Rainfall data not available prior to 1891")
  } else {
    if (any(yrs < 1960)) stop("Temperature data not available prior to 1960")
  }

  # Get date text used in file names - account for leap years
  mtxt<-ifelse(mnths<10,paste0("0",mnths),paste0("",mnths))
  daysofmonth<-c(31,28,31,30,31,30,31,31,30,31,30,31)
  mdays<-daysofmonth[mnths]
  mdays<-ifelse((yrs%%4==0 & mnths==2),29,mdays)

  # Create base url
  urlbase<-"https://dap.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.3.0.ceda/1km"

  # Create handle with header
  output<-c()
  for (v in varn){
    urls<-file.path(urlbase,v,"day","latest")
    fnames<-paste0(v,"_hadukgrid_uk_1km_day_",yrs,mtxt,"01-",yrs,mtxt,mdays,".nc")
    destfiles<-file.path(dir_out,fnames)
    dload_urls<-file.path(urls,fnames)
    output<-c(output,destfiles)
    for(n in 1:length(dload_urls)){
      print(paste('Downloading',v,'file',n))
      h <- new_handle(verbose = FALSE)
      handle_setheaders(h, .list=list("Authorization" = paste('Bearer',access_token)))
      curl_download(dload_urls[n],destfile=destfiles[n], handle = h)
    }
    #results<-multi_download(dload_urls,destfiles=destfiles,resume=TRUE,progress=TRUE)
  }
  # Checks downloaded files exist
  if(any(file.exists(output))==FALSE){
    warning('Download files not found!!!')
    print(output[which(file.exists(output)==FALSE)])
  }
  # Checks downloaded files >40MB
  info<-file.info(output)
  if(any(info$size<40000000)){
    warning('Downloaded files not expected size!!!')
    print(output[which(info$size<40000000)])
  }
  return(output)
}


#' Download UKCP18 climate data
#' @description
#' Using parameters the function will downloaded from ftp.ceda.ac.uk all available UCKP18 files containing relevant data to the user defined output directory
#' @param dir_out Output directory to which files are downloaded
#' @param startdate POSIXlt class defining start date of required timeseries
#' @param enddate POSIXlt class defining end date of required timeseries
#' @param collection text string defining UKCP18 collection, either 'land-gcm' or 'land-rcm'
#' @param domain text string defining UKCP18 domain, either 'uk' or 'eur'(land-rcm collection only) or 'global'
#' @param rcp text string iof RCP scenario to be downloaded
#' @param member vector of strings defining the climate model member to be downloaded. Available members vary between UKCP18 collections.
#' @param vars vector of strings corresponding to UKCP18 short variable names to download
#' @param download_dtm - if TRUE and orography (elevation) available will also download to dir_out
#' @param cedausr ceda username string
#' @param cedapwd ceda password sting
#' @param access_token ceda access token as provided by \href{https://services-beta.ceda.ac.uk/account/token/?_ga=2.221117848.1874942177.1727205775-405001952.1710779653}{Ceda access token generator}
#' @return downloads files to `dir_out`
#' @details start and end times are only used to identify decadal files containing data for this
#' time period - downloaded files are not cut to time period
#' @import curl
#' @export
#' @keywords download ukcp18
#' @examples
#'  \dontrun{
#' dir_ukcp <- tempdir()
#' cedausr<-"your_user_name"
#' cedapwd <- "your_password"
#' access_token<-"your_access_token"
#' download_ukcp18(dir_ukcp,as.POSIXlt('2018-05-01'),as.POSIXlt('2018-05-31'),'land-rcm','uk','rcp85',c('01'),c('tasmax','tasmin'),download_dtm=TRUE, cedausr,cedapwd,access_token)
#' }
download_ukcp18<-function(
    dir_out,
    startdate,
    enddate,
    collection=c('land-gcm','land-rcm'),
    domain=c('uk','eur','global'),
    rcp=c('rcp85', 'rcp26'),
    member=c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15',
             '16','17','18','19','20','21','22','23','24','25','26','27','28'),
    vars=c('clt','hurs','huss','pr','prsn','psl','rls','rss','tasmax','tasmin','uas','vas'),
    download_dtm=FALSE,
    cedausr,
    cedapwd,
    access_token
){
  # Check parameters
  if(!dir.exists(dir_out)) stop(paste("Output directory does not exist:",dir_out))
  collection<-match.arg(collection)
  domain<-match.arg(domain)
  rcp<-match.arg(rcp)
  member<-match.arg(member,several.ok=TRUE)
  vars<-match.arg(vars,several.ok=TRUE)
  if(collection=='land-rcm' & !any(member %in% c('01','04','05','06','07','08','09','10','11','12','13','15'))){
    warning('Invalid member for land-rcm - retricting to only those valid!!')
    member<-member[which(member %in% c('01','04','05','06','07','08','09','10','11','12','13','15'))]
  }
  if(collection=='land-gcm' & !domain %in% c('uk','global')){
    warning('Invalid area for global data - downloading global data!!')
    domain<-'global'
  }
  if(collection=='land-rcm' & !domain %in% c('uk','eur')){
    warning('Invalid area for regional data - downloading Europe data!!')
    domain<-'eur'
  }
  if(download_dtm & (collection!='land-rcm'|domain!='uk')){
    warning('DTM only downloadable land-rcm collection and uk domain - will not download!!! ')
    download_dtm<-FALSE
  }
  # Identify which decades are required
  decades<-.find_ukcp_decade(collection,startdate,enddate)

  # Get tiled data according to resolution
  if(collection=='land-rcm') res<-"12km"
  if(collection=='land-gcm') res<-"60km"

  url<-file.path("https://dap.ceda.ac.uk","badc","ukcp18","data",collection,domain,res,rcp)

  print('Parameters chosen:')
  print(collection)
  print(domain)
  print(rcp)
  print(member)
  print(vars)
  print(decades)
  print(paste('Downloading from',url))

  # Download files - loop over model runs and vars - use lapply to download all decades
  output<-c()
  for(run in member){
    for(v in vars){
      fnames<-paste0(v,'_',rcp,'_',collection,'_',domain,'_',res,'_',run,'_day_',decades,'.nc')
      dload_urls<-file.path(url,run,v,'day','latest',fnames)
      destfiles<-file.path(dir_out,fnames)
      output<-c(output,destfiles)
      print(length(destfiles)==length(dload_urls))
      for(n in 1:length(dload_urls)){
        print(paste('Downloading',run,v,'file',n))
        h <- new_handle(verbose = FALSE)
        handle_setheaders(h, .list=list("Authorization" = paste('Bearer',access_token)))
        curl_download(dload_urls[n],destfile=destfiles[n], handle = h)
      }
      #h <- curl::new_handle()
      #curl::handle_setopt(h, userpwd = paste0(cedausr,":",cedapwd))
      #success_df<-curl::multi_download(urls=dload_urls, destfiles=destfiles, userpwd = paste0(cedausr,":",cedapwd))
      #if(any(success_df$success==FALSE)) print(paste('UNSUCCESSFUL file downloads:',fnames[which(!success_df$success)])) else print('All downloads successful')
      #if(nrow(fullreport_df)==0) fullreport_df<-success_df else fullreport_df<-rbind(fullreport_df,success_df)
    }
  }
  # Download orography if requested
  if(download_dtm){
    print('Downloading orography for UK rcm...')
    fname<-'orog_land-rcm_uk_12km_osgb.nc'
    orogfile<-file.path(dir_out,fname)
    orog_url<-paste0("https://dap.ceda.ac.uk/badc/ukcp18/data/land-rcm/ancil/orog/",fname)
    h <- new_handle(verbose = FALSE)
    handle_setheaders(h, .list=list("Authorization" = paste('Bearer',access_token)))
    curl_download(orog_url,destfile=orogfile, handle = h)
  }
  # Checks expected downloaded files exist
  if(any(file.exists(output))==FALSE){
    warning('Download files not found!!!')
    print(output[which(file.exists(output)==FALSE)])
  }
  # Checks downloaded files >100MB
  info<-file.info(output)
  if(any(info$size<100000000)){
    warning('Downloaded files not expected size!!!')
    print(output[which(info$size<100000000)])
  }
  if(download_dtm) output<-c(destfiles,orogfile)
  return(output)
}


#' @title Download UKCP18rcm driven sea surface temperature predictions
#' @description For description of the modelling used to produce sea surface temperatures then see:
#'  `https://www.metoffice.gov.uk/binaries/content/assets/metofficegovuk/pdf/research/ukcp/ukcp18-guidance-data-availability-access-and-formats.pdf`
#' @param dir_out - local directory to which files written
#' @param startdate - start date as POSIXlt
#' @param enddate  - end date as POSIXlt
#' @param modelruns - vector of two character strings of model numbers
#' @param cedausr ceda username string
#' @param cedapwd ceda password sting
#' @param access_token ceda access token as provided by \href{https://services-beta.ceda.ac.uk/account/token/?_ga=2.221117848.1874942177.1727205775-405001952.1710779653}{Ceda access token generator}
#'
#' @return downloaded file paths
#' @export
#' @keywords download ukcp18
#' @examples
#' \dontrun{
#'  cedausr<-"your_user_name"
#'  cedapwd <- "your_password"
#'  access_token<-"your_access_token"
#'  startdate<-as.POSIXlt('2017/12/31')
#'  enddate<-as.POSIXlt('2018/12/31')
#'  modelruns<-c('01')
#'  dir_out<-tempdir()
#'  download_ukcpsst(dir_out,startdate,enddate,modelruns, cedausr,cedapwd,access_token)
#'  }
download_ukcpsst<-function(
    dir_out,
    startdate,
    enddate,
    modelruns=c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15',
                '16','17','18','19','20','21','22','23','24','25','26','27','28'),
    cedausr,
    cedapwd,
    access_token
){
  # Check parameters
  modelruns<-match.arg(modelruns,several.ok=TRUE)
  # Get member ID used in file names from mesoclim lookup table
  memberid<-ukcp18lookup$PP_ID[which(ukcp18lookup$Member_ID %in% modelruns)]
  if("" %in% memberid) warning(paste("Model members",modelruns[which(memberid=="")],"NOT available for sea surface temperature - ignoring!!"))
  memberid<-memberid[which(memberid!="")]

  if(!dir.exists(dir_out))(stop(paste("Directory",dir_out,"does not exist!!!")))
  if(class(startdate)[1]!="POSIXlt" | class(enddate)[1]!="POSIXlt") stop("Date parameters NOT POSIXlt class!!")

  # Derive months of data required - will output month before and after start and end dates for interpolation
  start<-startdate %m-% months(1)
  end<-enddate %m+% months(1)
  yrs<-unique(c(year(start):year(end)))

  output<-c()
  for(id in memberid){
    url<-paste0("https://dap.ceda.ac.uk/badc/deposited2023/marine-nwsclim/NWSPPE/",id,"/annual")
    fnames<-paste0("NWSClim_NWSPPE_",id,"_",yrs,"_gridT.nc")
    dload_urls<-file.path(url,fnames)
    destfiles<-file.path(dir_out,fnames)
    output<-c(output,destfiles)
    for(n in 1:length(dload_urls)){
      print(paste('Downloading',id,'file',n))
      h <- new_handle(verbose = FALSE)
      handle_setheaders(h, .list=list("Authorization" = paste('Bearer',access_token)))
      curl_download(dload_urls[n],destfile=destfiles[n], handle = h)
    }
  }
  # Checks expected downloaded files exist
  if(any(file.exists(output))==FALSE){
    warning('Download files not found!!!')
    print(output[which(file.exists(output)==FALSE)])
  }
  # Checks downloaded files >40MB
  info<-file.info(output)
  if(any(info$size<40000000)){
    warning('Downloaded files not expected size!!!')
    print(output[which(info$size<40000000)])
  }
  return(output)
}

#' @title Download 1km albedo data
#' @description The function `download_globalbedo` downloads monthly GlobAlbedo tiled
#' albedo data from Jasmin available for the years 1998-2011. Whole of the UK covered by default tiles. For further information see:
#' \url{http://www.globalbedo.org/index.php}
#' @param dir_out directory to which to save data
#' @param year numeric vector indicating the years for which data are required
#' @param month string vector  indicating the month for which data are required
#' @param tiles vector of strings of the file name in the format 'h##v##'. Whole of the UK covered by default tiles.
#' @details No username or password currently required.
#' @return database detailing download details and success of all requested files
#' @import curl
#' @export
#' @keywords download
download_globalbedo<-function(dir_out,
                              tiles=c('h17v03','h17v04','h18v03'),
                              months=c('01','02','03','04','05','06','07','08','09','10','11','12'),
                              years=seq(1998,2011),
                              url='http://gws-access.jasmin.ac.uk'
){
  # Check parameters
  if(!dir.exists(dir_out)) stop(paste("Output directory does not exist:",dir_out))
  months<-match.arg(months,several.ok=TRUE)
  years<-match.arg(years,several.ok=TRUE)
  tiles<-match.arg(tiles,several.ok=TRUE)

  path<-'public/qa4ecv/albedo/netcdf_cf/1km/monthly'
  fullreport_df<-data.frame()

  for(t in tiles){
    for(y in years){
      print(paste('Downloading files for ',t,'-',y))
      fnames<-paste0('GlobAlbedo.merge.albedo.',y,months,'.',t,'.nc')
      dload_urls<-file.path(url,path,t,y,fnames)
      destfiles<-file.path(dir_out,fnames)

      success_df<-curl::multi_download(urls=dload_urls, destfiles=destfiles )
      if(any(success_df$success==FALSE)) print(paste('UNSUCCESSFUL file downloads:',fnames[which(!success_df$success)])) else print('All downloads successful')
      if(nrow(fullreport_df)==0) fullreport_df<-success_df else fullreport_df<-rbind(fullreport_df,success_df)
    }
  }
  return(fullreport_df)
}
# ------------------------------------ ------------------------------------ ------------------
# ------------------------------------ ERA5 DOWNLOADS ------------------------------------
# ------------------------------------ ------------------------------------ ------------------
#' Download ERA5 reanalysis data
#'
#' @param dir_out - Local output directory
#' @param file_out - Output filename
#' @param year - year of data to retrieve
#' @param area - bounding box with order c(ymax,xmin,ymin,xmax)
#' @param variables - Long names of variables to retrieve
#' @param era5_user - cds user name/email
#' @param era5_key - cds api key
#'
#' @return filepath of output
#' @details Downloads a whole single year of reanalysis-era5-single-levels as netcdf file.
#' Variables chosen and size of area defaults match UK/Eire and variables required for microclimate modelling.
#' Limits on download file may restrict size of area and/or number of variables that can be downloaded to single file.
#' @import curl
#' @import ecmwfr
#' @export
#'
#' @examples
#' \dontrun{
#' dir_out <- tempdir()
#' cds_user<-"your_user_name"
#' cds_key<-"your_key_string"
#' download_era5(dir_out,file_out='test.nc',year=2020,area = c(55,-3,49,3),variables =c("2m temperature"),era5_user=cds_user,era5_key=cds_key)
#' }
download_era5<-function(dir_out,
                        file_out,
                        year,
                        area = c(61, -11.25, 48.75, 2.5),
                        variables = c("10m u-component of wind",
                                      "10m v-component of wind",
                                      "2m dewpoint temperature",
                                      "2m temperature",
                                      "Mean sea level pressure",
                                      "Surface pressure",
                                      "Total precipitation",

                                      "Surface solar radiation downwards",
                                      "Total sky direct solar radiation at surface",

                                      "Mean surface downward short-wave radiation flux",
                                      "Mean surface direct short-wave radiation flux",
                                      "Mean surface downward long-wave radiation flux",
                                      "Mean surface net long-wave radiation flux" ),
                        era5_user,
                        era5_key
){

  wf_set_key(key=era5_key,user= era5_user)

  request <- list(
    dataset_short_name = "reanalysis-era5-single-levels",
    product_type = "reanalysis",
    variable = variables,
    year = year,
    month =c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'),
    day = c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12',
            '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24',
            '25', '26', '27', '28', '29', '30', '31'),
    time =c('00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00', '12:00', '13:00', '14:00', '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00', '21:00', '22:00', '23:00'),
    data_format = "netcdf",
    download_format = "unarchived",
    area = area,
    target = file_out
  )
  file<-wf_request(request = request,
                   user = era5_user,
                   transfer = TRUE,
                   path = dir_out)
  return(file)
}

#' Download ancillary era5 data (geopotential and lansea mask)
#'
#' @param dir_out - directory to download files to
#' @param area - lat lon of area
#' @param variables one or both of geopotential and land_sea_mask
#' @param era5_user CDS user name
#' @param era5_key CDS API token
#'
#' @return filepaths of outputs
#' @details Downloads a single instance of ancillary dat for reanalysis-era5-single-levels as netcdf files.
#' @import curl
#' @import ecmwfr
#' @export
#'
#' @examples
#' \dontrun{
#' dir_out <- tempdir()
#' cds_user<-"your_user_name"
#' cds_key<-"your_key_string"
#' download_ancillary_era5(dir_out,area = c(55,-3,49,3),variables =c("geopotential","land_sea_mask" ),era5_user=cds_user,era5_key=cds_key)
#' }
download_ancillary_era5<-function(dir_out,
                                  area = c(61, -11.25, 48.75, 2.5),
                                  variables = c("geopotential",
                                                "land_sea_mask" ),
                                  era5_user,
                                  era5_key
){

  wf_set_key(key=era5_key,user= era5_user)
  for(v in variables){
    request <- list(
      dataset_short_name = "reanalysis-era5-single-levels",
      product_type = "reanalysis",
      variable = v,
      year = 2020,
      month =c('01'),
      day = c('01'),
      time =c('00:00'),
      data_format = "netcdf",
      download_format = "unarchived",
      area = area,
      target = paste0("era5_",v,".nc")
    )
    file<-wf_request(request = request,
                     user = era5_user,
                     transfer = TRUE,
                     path = dir_out)
  }
  return(file)
}
