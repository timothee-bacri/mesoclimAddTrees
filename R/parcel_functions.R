#' @title Calculates climate timeseries for vector parcels
#' @description Calculates timeseries of climate variables for each polygon in parcels using weighted means of overlapping grid cells.
#'
#' @param climdata - list of SpatRasters of downscaled mesoclimate variables as output by `spatialdownscale`
#' @param parcels - sf or vect describing polygons for which climate data is required
#' @param id - string of variable name in parcels corresponding to parcel identification field
#'
#' @return nested list of climate variable timeseries by parcel.
#' @export
#' @import terra
#' @keywords postprocess data
create_parcel_list<-function(climdata,parcels,id='gid',
                             input_names=c("tmax", "tmin","swrad","lwrad","relhum","pres","prec", "windspeed"),
                             output_names=c("tmax", "tmin","swdown","lwdown","relhum","pres","prec", "windspeed"),
                             roundings=c(2,2,1,1,1,1,2,3)){
  if(any(!input_names %in% names(climdata))) stop('Input name NOT found in climate dataset provided')
  if(length(input_names)!=length(output_names)) stop('Different number of input and output names!!!')
  if(length(roundings)!=length(output_names)) stop('Different number of rounding values and output names!!!')
  dtmf<-climdata$dtm

  # Create empty df for a parcel
  rws<-length(climdata$tme)
  parcel_df<-as.data.frame(matrix(ncol=1+length(output_names), nrow=rws))
  names(parcel_df)<-c("date",output_names)
  parcel_df$date=as.Date(as.character(climdata$tme))

  # Create a list of parcel_dfs
  parcel_list<- lapply(1:length(parcels), function(x) parcel_df)
  names(parcel_list)<-values(parcels)[,id]

  # Extract weighted means for each parcel by climate variable
  for (n in 1:length(input_names)){
    vin<-input_names[n]
    vout<-output_names[n]
    rnd<-roundings[n]
    print(vin)
    r<-climdata[[vin]]
    vals<-t(terra::extract(r,parcels,fun=mean,weights=TRUE,raw=TRUE,na.rm=TRUE, ID=FALSE)) # matrix of timestep x parcel
    for(n in 1:length(parcels)) parcel_list[[n]][[vout]]<-round(vals[,n],rnd)
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
