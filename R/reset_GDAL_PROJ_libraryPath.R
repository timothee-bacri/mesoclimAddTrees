# Reset gdal library for R
# see: https://github.com/rspatial/terra/issues/1378

Sys.getenv("PROJ_LIB")
# plib = "C:\\OSGeo4W64\\share\\proj"

prj <- system.file("proj", package = "sf")[1]
Sys.setenv("PROJ_LIB" = prj)
Sys.getenv("PROJ_LIB")

# THEN LOAD LIBRARIES eg terra, sf

# and perhaps set it back when done so that you can use Postgres
# Sys.setenv("PROJ_LIB" = "C:\\OSGeo4W64\\share\\proj")
library(devtools)
load_all()

cedausr<-'jrmosedale'
cedapwd<-"T${5Q9Z<9ef'"
