#' @title Applies bias correction to UK era5 temperature data
#' @description The function `correct_era5temps` applies automatic bias correction to era5
#' temperature data data to correct for unaturally low diurnal temperature
#' fluctuations in coastal grid cells.
#' @param era5hourly a stacked SpatRast of hourly ERA5 temperature data for any part
#' of the UK.
#' @param era5correctmodels a list of model correction coefficients for each UK
#' grid cell. Available: https://universityofexeteruk-my.sharepoint.com/:f:/g/personal/i_m_d_maclean_exeter_ac_uk/EjxJdJq0MClFqUU3LATlcnEBo1cGFiUAxqLQALNNxvdZaw?e=wLR2Rf
#' Zenodo link to be added dreckly.
#' @return a SpatRast of bias corrected temperature data
#' @details `era5correctmodels` was derived by applying [biascorrect()] to a 2018
#' dataset of era5 temperature data, calibrating against Met office data
#' @import terra
#' @import mgcv
#' @import mesoclim
#' @export
#' @keywords biascorrect
#' @rdname correct_era5temps
correct_era5temps<-function(era5hourly,era5correctmodels) {
  # Check whether in Kelvin or degrees
  me<-mean(.is(era5hourly[[1]]),na.rm=TRUE)
  if (me > 150) era5hourly<-era5hourly-273.15
  # Calculate tmx, tmn and dtr
  tmx_era5<-hourtodayCpp(.is(era5hourly), "max")
  tmn_era5<-hourtodayCpp(.is(era5hourly), "min")
  dtr_era5<-tmx_era5-tmn_era5
  # Match cells in era5hourly to model list index
  # get x & ys of model list index
  rindex<-rast(era5correctmodels$indexr)
  xy1<-data.frame(xyFromCell(rindex, 1:ncell(rindex)))
  z<-as.vector(extract(rindex,xy1)[,2])
  xyz1 <- cbind(xy1,z)
  s<-which(is.na(z)==F)
  xyz1<-xyz1[s,]
  xyz1$x<-round(xyz1$x*4,0)/4
  xyz1$y<-round(xyz1$y*4,0)/4
  # get x & ys of model era5hourly grid cells
  xy2<-data.frame(xyFromCell(era5hourly[[1]], 1:ncell(era5hourly[[1]])))
  v<-c(1:(dim(xy2)[1]))
  xyz2 <- cbind(xy2,v)
  xyz2$x<-round(xyz2$x*4,0)/4
  xyz2$y<-round(xyz2$y*4,0)/4
  vindex<-rast(xyz2)
  mindex<-.is(vindex)
  # merge datasets
  xy<-merge(xyz1,xyz2,by=c("x","y"),all=TRUE)
  adtr<-array(NA,dim=dim(dtr_era5))
  atmn<-adtr
  for (i in 1:dim(adtr)[1]) {
    for (j in 1:dim(adtr)[2]) {
      v1<-dtr_era5[i,j,]
      v2<-tmn_era5[i,j,]
      # get model
      v<-mindex[i,j]
      s<-which(xy$v==v)
      index<-xy$z[s]
      if (is.na(index) == FALSE) {
        mod1<-era5correctmodels$dtrmod[[index]]
        mod2<-era5correctmodels$tmnmod[[index]]
        adtr[i,j,] <- predict.gam(mod1, newdata = data.frame(v2 = v1))
        atmn[i,j,] <- predict.gam(mod2, newdata = data.frame(v2 = v2))
      }
    }
  }
  # Apply dtr and min correction to original datasets
  # ~~ Calculate hourly dtr fractions
  dtr_era5h<-.ehr(dtr_era5)
  tmn_era5h<-.ehr(tmn_era5)
  tc<-.is(era5hourly)
  tfrac<-(tc-tmn_era5h)/dtr_era5h
  # Apply correct data
  adtrh<-.ehr(adtr)
  atmnh<-.ehr(atmn)
  tcn<-(tfrac*adtrh)+atmnh
  tcn<-.rast(tcn,era5hourly)
  return(tcn)
}

#' @title Blends Met Office and ERA5 data to produce hourly 1km resolution temperature data
#' @description The function `blendtemp_hadukera5` ERA5 data to 1 km grid resoltuion,
#' calculates the diurnal cycle in each grid cell, and then adjusts this by the
#' maximum and minimum daily temperatures in the one km met office data.
#' @param tasmin a stacked SpatRaster of haduk daily minimum temperatures (deg C)
#' @param tasmax a stacked SpatRaster of haduk daily maximum temperatures (deg C)
#' @param era5t2m a stacked SpatRaster of hourly ERA5 temperatures (deg C or K)
#' @import mesoclim
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @export
#' @keywords temporal preprocess
blendtemp_hadukera5<-function(tasmin,tasmax,era5t2m) {
  d1<-dim(tasmin)
  d2<-dim(tasmax)
  d3<-dim(era5t2m)
  if (sum(d1-d2) != 0) stop("dims of tasmin and tasmac must match")
  if (d3[3] != (d2[3]*24)) stop("Hours in era5 must match days in tasmin")
  # Reproject era5
  era5t2m<-project(era5t2m,tasmin)
  # met office dtr
  dtr<-tasmax-tasmin
  # era5 dtr and min
  era5min<-.hourtodayCpp(.is(era5t2m),"min")
  era5max<-.hourtodayCpp(.is(era5t2m),"max")
  era5dtr<-era5max-era5min
  era5minh<-.ehr(era5min)
  era5dtrh<-.ehr(era5dtr)
  era5frac<-(.is(era5t2m)-era5minh)/era5dtrh
  # met office hourly
  tasminh<-.ehr(.is(tasmin))
  dtrh<-.ehr(.is(dtr))
  tch<-.rast(xx<-era5frac*dtrh+tasminh,tasmin)
  return(tch)
}

#' @title Applies biascorrectukcpone (precipcorrect) to all ukcp variables for one tile and model run
#' and saves data to disk as nc files per decade
#' @param pathtoera - directory with 2018 daily era5 data (as returned by era5todaily)
#' @param pathtoukcp18 - directory with 2018 ukcp data (as returned by crop2018UKCP)
#' @param pathtoukcpdecade - directory with decadal ukcp data (as returned by cropandsortUKCPdecade)
#' @param pathout - directory in which to save corrected data
#' @param decade - the decade for which data are required (1 for 2010-2019, 2 for 2020-2019 etc)
#' @param modelrun - a numeric value (1 is convertyed to 01) indicating the model run. Used for reading in
#' data, so file naming and foldr convention assumed to match that of data supplied via dropbox
#' saves bias corrected daily ukcp data as compressed ncdf4 files
#' @import mesoclim
#' @noRd
biascorrectukcpall<-function(pathtoera,pathtoukcp18,pathtoukcpdecade,pathout,decade,modelrun) {
  # file names
  mtxt<-ifelse(modelrun<10,paste0("0",modelrun),paste0("",modelrun))
  to<-paste0(2000+decade*10,"_",2000+decade*10+9)
  byr<-2000+decade*10
  # Temperature
  # ** era5
  fi<-paste0(pathtoera,"tasmax.tif")
  tmx<-rast(fi)
  fi<-paste0(pathtoera,"tasmin.tif")
  tmn<-rast(fi)
  edtr<-tmx-tmn
  etme<-(tmx+tmn)/2
  # ** ukcp fit
  fi<-paste0(pathtoukcp18,"tasmax_2018.tif")
  tmx<-rast(fi)
  fi<-paste0(pathtoukcp18,"tasmin_2018.tif")
  tmn<-rast(fi)
  udtr<-tmx-tmn
  utme<-(tmx+tmn)/2
  # ** ukcp apply
  fi<-paste0(pathtoukcpdecade,"tasmax_",mtxt,"_",to,".tif")
  tmx<-rast(fi)
  fi<-paste0(pathtoukcpdecade,"tasmin_",mtxt,"_",to,".tif")
  tmn<-rast(fi)
  uadtr<-tmx-tmn
  uatme<-(tmx+tmn)/2
  dtrc<-biascorrectukcpone(edtr, udtr, uadtr, rangelims = 1.1)
  tmec<-biascorrectukcpone(etme, utme, uatme, rangelims = 1.1)
  tmxc<-tmec+0.5*dtrc
  tmnc<-tmec-0.5*dtrc
  fo1<-paste0(pathout,"tasmax_",mtxt,"_",to,".nc")
  fo2<-paste0(pathout,"tasmin_",mtxt,"_",to,".nc")
  savenc(tmxc,byr,"tasmax","Maximum air temperature at 2 m","deg C",fo1)
  savenc(tmnc,byr,"tasmin","Minimum air temperature at 2 m","deg C",fo2)
  # Pressure
  fi<-paste0(pathtoera,"psl.tif")
  r1<-rast(fi)
  fi<-paste0(pathtoukcp18,"psl_2018.tif")
  r2<-rast(fi)
  fi<-paste0(pathtoukcpdecade,"psl_",mtxt,"_",to,".tif")
  r3<-rast(fi)
  psl<-biascorrectukcpone(r1, r2, r3, rangelims = 1.1)
  fo<-paste0(pathout,"psl_",mtxt,"_",to,".nc")
  savenc(psl,byr,"psl","Surface level pressure","kPa",fo)
  # Humidity
  fi<-paste0(pathtoera,"huss.tif")
  r1<-rast(fi)
  fi<-paste0(pathtoukcp18,"huss_2018.tif")
  r2<-rast(fi)
  fi<-paste0(pathtoukcpdecade,"huss_",mtxt,"_",to,".tif")
  r3<-rast(fi)
  huss<-biascorrectukcpone(r1, r2, r3, rangelims = 1.1)
  h<-as.array(huss)
  relh<-suppressWarnings(converthumidity(h, intype = "specific", outtype = "relative",
                                         tc = as.array(tmec), pk = as.array(psl)))
  relh[relh>100]<-100
  relh[relh<20]<-20
  relh<-rast(relh)
  ext(relh)<-ext(huss)
  crs(relh)<-crs(huss)
  fo<-paste0(pathout,"relhum_",mtxt,"_",to,".nc")
  savenc(relh,byr,"relhum","Relative humidity at 2 m","Percentage",fo)
  # Wind speed and direction
  # ** u wind
  ufi<-paste0(pathtoera,"uas.tif")
  ur1<-rast(ufi)
  ufi<-paste0(pathtoukcp18,"uas_2018.tif")
  ur2<-rast(ufi)
  ufi<-paste0(pathtoukcpdecade,"uas_",mtxt,"_",to,".tif")
  ur3<-rast(ufi)
  # ** v wind
  vfi<-paste0(pathtoera,"vas.tif")
  vr1<-rast(vfi)
  vfi<-paste0(pathtoukcp18,"vas_2018.tif")
  vr2<-rast(vfi)
  vfi<-paste0(pathtoukcpdecade,"vas_",mtxt,"_",to,".tif")
  vr3<-rast(vfi)
  # ** Calculate wind speed
  ws1<-sqrt(ur1^2+vr1^2)
  ws2<-sqrt(ur2^2+vr2^2)
  ws3<-sqrt(ur3^2+vr3^2)
  wsc<-biascorrectukcpone(ws1, ws2, ws3, rangelims = 1.1)
  # ** Adjust to 2 m
  wsc<-wsc*(4.87/log(67.8*10-5.42))
  wsc[wsc<0.001]<-0.001
  wdir<-(atan2(ur3,vr3)*(180/pi)+180)%%360
  # ** Save files
  fo<-paste0(pathout,"winds_",mtxt,"_",to,".nc")
  savenc(wsc,byr,"winds","Wind speed at 2 m","m/s",fo)
  fo<-paste0(pathout,"windd_",mtxt,"_",to,".nc")
  savenc(wdir,byr,"windd","Wind direction","degrees",fo)
  # ** Sky emmisivity
  fi<-paste0(pathtoera,"skyem.tif")
  r1<-rast(fi)
  fi<-paste0(pathtoukcp18,"skyem_2018.tif")
  r2<-rast(fi)
  fi<-paste0(pathtoukcpdecade,"skyem_",mtxt,"_",to,".tif")
  r3<-rast(fi)
  skyem<-biascorrectukcpone(r1, r2, r3, rangelims = 1.1)
  skyem[skyem>1]<-1
  skyem[skyem<0]<-0
  fo<-paste0(pathout,"skyem_",mtxt,"_",to,".nc")
  savenc(skyem,byr,"skyem","Sky emissivity","0-1",fo)
  # ** Shortwave radiation
  fi<-paste0(pathtoera,"rss.tif")
  r1<-rast(fi)
  fi<-paste0(pathtoukcp18,"rss_2018.tif")
  r2<-rast(fi)
  fi<-paste0(pathtoukcpdecade,"rss_",mtxt,"_",to,".tif")
  r3<-rast(fi)
  rss<-biascorrectukcpone(r1, r2, r3, rangelims = 1.1)
  rss[rss<0]<-0
  fo<-paste0(pathout,"rss_",mtxt,"_",to,".nc")
  savenc(rss,byr,"rss","Total downward shortwave radiation","W/m**2",fo)
  # ** Precipitation
  fi<-paste0(pathtoera,"pr.tif")
  r1<-rast(fi)
  fi<-paste0(pathtoukcp18,"pr_2018.tif")
  r2<-rast(fi)
  fi<-paste0(pathtoukcpdecade,"pr_",mtxt,"_",to,".tif")
  r3<-rast(fi)
  pr<-precipcorrect(r1,r2,r3)
  fo<-paste0(pathout,"pr_",mtxt,"_",to,".nc")
  savenc(pr,byr,"pr","Precipitation rate","mm/day",fo)
}

