#' @title Calculates climate timeseries for vector parcels
#' @description Calculates timeseries of climate variables for each polygon in parcels using weighted means of overlapping grid cells.
#'
#' @param climdata - list of SpatRasters of downscaled mesoclimate variables as output by `spatialdownscale`
#' @param parcels - sf or vect describing polygons for which climate data is required
#' @param id - string of variable name in parcels corresponding to parcel identification field
#'
#' @return list of dataframes of climate variable timeseries by parcel.
#' @export
#' @importFrom sf st_as_sf
#' @importFrom exactextractr exact_extract
#' @keywords postprocess data
create_parcel_list<-function(climdata,parcels,id='gid',
                             input_names=c("tmax", "tmin","swrad","lwrad","relhum","pres","prec", "windspeed"),
                             output_names=c("tmax", "tmin","swdown","lwdown","relhum","pres","prec", "windspeed"),
                             roundings=c(2,2,1,1,1,1,2,3)){
  if(any(!input_names %in% names(climdata))) stop('Input name NOT found in climate dataset provided')
  if(length(input_names)!=length(output_names)) stop('Different number of input and output names!!!')
  if(length(roundings)!=length(output_names)) stop('Different number of rounding values and output names!!!')
  parcels<-st_as_sf(parcels)
  dtmf<-climdata$dtm

  # Create empty df for a parcel
  rws<-length(climdata$tme)
  parcel_df<-as.data.frame(matrix(ncol=1+length(output_names), nrow=rws))
  names(parcel_df)<-c("date",output_names)
  parcel_df$date=as.Date(as.character(climdata$tme))

  # Create a list of parcel_dfs
  parcel_list<- lapply(1:nrow(parcels), function(x) parcel_df)
  names(parcel_list)<-parcels[[id]]

  # Extract weighted means for each parcel by climate variable
  for (n in 1:length(input_names)){
    vin<-input_names[n]
    vout<-output_names[n]
    rnd<-roundings[n]
    print(vin)
    r<-climdata[[vin]]
    #vals<-t(terra::extract(r,parcels,fun=mean,weights=TRUE,raw=TRUE,na.rm=TRUE, ID=FALSE)) # matrix of timestep x parcel
    vals<-round(exact_extract(r,parcels,'mean', weights='area'),rnd)
    #for(n in 1:length(parcels)) parcel_list[[n]][[vout]]<-round(vals[,n],rnd)
    for(n in 1:length(parcels)) parcel_list[[n]][[vout]]<-as.numeric(vals[n,])
  }

  return(parcel_list)
}


#' @title Write parcel climate .csv files
#' @description Writes one file per parcel of climate timeseries with option to overwrite or append to existing files.
#'
#' @param parcel_list - list of parcel names and climate variables as output by `create_parcel_list`
#' @param dir_out - directory to write to
#' @param overwrite - string describing how to respond to existing files where 'none' stops if any file exists,
#' 'overwrite' replaces existing files and 'append' adds results to existing files.
#'
#' @return NA - writes .csv files
#' @importFrom utils write.table
#' @export
#' @keywords postprocess data
write_parcels<-function(parcel_list, dir_out, overwrite=c('none','append','replace')){
  overwrite<-match.arg(overwrite)
  parcel_names<-names(parcel_list)
  # If overwrite= none then check if any output files already exist and stop if true
  if(overwrite=='none'){
    fnames<-file.path( dir_out, paste0("parcel_",parcel_names[n],"_mesoclim.csv") )
    if(any(file.exists(fnames))) stop('Output file exists already!!! Use overwrite=replace to replace.')
  }
  # Write files
  for(n in 1:length(parcel_list)){
    fout<-file.path( dir_out, paste0("parcel_",parcel_names[n],"_mesoclim.csv") )
    paste(fout)
    if(overwrite=='replace') write.table(parcel_list[[n]], fout, sep = ",", row.names=FALSE, col.names = TRUE)
    if(overwrite=='append') write.table(parcel_list[[n]], fout,sep = ",", row.names=FALSE, col.names = FALSE, append = T)
  }
}

#' Get parcel variable
#' @description - get parcel values of chosen climate variable as mean/min/max values over all timesteps.
#' Output suitable for easy mapping using `plot_parcel_var` function
#' @param climdata - named list of climate timeseriess
#' @param var - name of climdata element to be got
#' @param parcels - vect or sf dataframe of parcel polygons
#' @param id - unique id variable in parcels
#' @param stat - summary statistic (min,max or mean) to be used to return single raster layer
#'
#' @return sf dataframe of parcels and variable values - suitable for mapping
#' @export
#' @importFrom terra app
#' @importFrom sf st_as_sf
#' @importFrom exactextractr exact_extract
#' @examples
get_parcel_var<-function(climdata,var, parcels,id='gid', stat=c('mean','min','max') ){
  if(!var %in% names(climdata)) stop("Variable not found in climdata provided!!!")
  if (class(climdata[[var]])[1]!="SpatRaster") stop("Date is not a spatraster!!!")
  results_sf<-st_as_sf(parcels)

  # Calculate stat across all layers
  stat_r<-app(climdata[[var]],fun=stat)

  # Extract variable
  results_sf[[var]]<-exact_extract(stat_r,results_sf,'mean', weights='area')
  return(results_sf)
}

#' Plot parcel variable as leaflet map
#'
#' @param parcels_sf - vect or sf object of parcels
#' @param plotvar - name of variable column to be plotted in parcels_sf
#' @param idvar - id variable in parcels_sf
#'
#' @return leaflet map object
#' @export
#' @import leaflet
#' @importFrom sf st_transform st_coordinates st_centroid st_union st_make_valid
#' @examples
plot_parcel_var<-function(parcels_sf, plotvar='tmax', idvar='gid'){
  minval<-floor(min(parcels_sf[[plotvar]]))
  maxval<-ceiling(max(parcels_sf[[plotvar]]))

  # Palette
  pal_brks <- seq(minval, maxval, by=(maxval-minval)/9)
  pal<-colorBin('BuGn',domain=c(pal_brks[1],pal_brks[length(pal_brks)+1]),bins=pal_brks,na.color="transparent",reverse=FALSE)

  # Popup info
  popup_text <- paste("<strong>Parcel ID: </strong>",  parcels_sf[[idvar]],
                       "Temp =",round(parcels_sf[[plotvar]],1))

  # Map
  parcels_sf<-st_transform(parcels_sf,4326)
  centre_pt<-st_coordinates(st_centroid(st_union(st_make_valid(parcels_sf))))

  varmap<-leaflet(options=list(minZoom=8)) %>%
    setView(lng = centre_pt[1], lat = centre_pt[2], zoom = 12) %>%
    addProviderTiles(providers$OpenStreetMap.Mapnik,group="Streetmap") %>%
    addPolygons(data=parcels_sf$geometry,col='black',fillColor=pal(parcels_sf[[plotvar]]),
                weight = 0.5, fillOpacity = 1, opacity = 1,
                popup = popup_text,
                highlightOptions = highlightOptions(color = "yellow",weight = 2,
                                                    fillColor = 'yellow', bringToFront = TRUE)) %>%
  addLegend(position = "bottomleft", pal=pal,values=c(pal_brks[1],pal_brks[length(pal_brks)+1]))
  return(varmap)
}

