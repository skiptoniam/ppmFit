#'@title Function for generating a spatial thinned point process to reflect effort.
#'@name effortOffset
#'@rdname effortOffset
#'@description This function helps you set up a spatial thinning process (the probability that a cell has been surveyed) from known
#'presence-only occurrence data and a defined spatial window. The approach returns a thinning probability which is equad to: pi = 1 - exp(-lambda)
#'where lambda is the observed number of counts per cell. Clearly the cell resolution will effect the lambda values, so
#'make sure the window is the same spatial scale and resolution as the model you intend to fit.
#'@param presences A two column data.frame or matrix of presence locations. Typically we might used the presences for all species that fall in a particular taxa (say all reptiles).
#'@param window A spatial window at the same resolution and extent as the data is to be modelled.
#'@param thin.prob returned the plug in thinned probability for to use in the offset.
#'@param buffer A spatial distance (on the scale of the raster resolution) to buffer around each point. This will count the number of points per cell with the size of the point buffered to a radius of 'buffer' size.
#'@param eps A small value to use as the zero for the thinned probability.
#'@author Skipton Woolley
#'@references Cressie, N. 1992. Statistics for spatial data. John Wiley & Sons..
#'@importFrom terra rasterize rasterizeWin mask app
#'@export
#'@examples
#'library(ppmData)
#'path <- system.file("extdata", package = "ppmData")
#'lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#'win <- rast(lst[1])
#'presences <- pres <- as.matrix(snails[,1:2])
#'thpp <- effortOffset(pres,win)

## test on big data set
# library(terra)
# library(sp)
# df <- read.csv("dev/reptilia_2022-10-06.csv")
# ausProj <- "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
# wgs84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
# xy84 <- sp::SpatialPoints(df[,c("longitude","latitude")], proj4string=CRS(wgs84))
# x <- spTransform(xy84, ausProj)
# presences <- coordinates(x)
# window <- rast("dev/test.tif")
# t1 <- effortOffset(presences,window)
# t2 <- effortOffset(presences,window,buffer=2500)

effortOffset <- function(presences, window, thin.prob=TRUE, buffer=NULL, eps=1e-9){

  if(missing(presences) | !is.matrix(presences)){
    stop("please make sure that presences is included as a two column matrix.")
  }

  if(missing(window)){
    stop("please make sure window is included. Window should be a terra spatial raster.")
  }

  if(is.null(buffer)){
    rlambda <- terra::rasterize(presences, window, fun="count")
  } else {
    message("'buffer' should be a buffer distance that is appropriate for the resolution of the raster.")
    rlambda <- terra::rasterizeWin(data.frame(presences,z=1), window, fun="count", pars=buffer, win="circle")
  }

  if(thin.prob){
    rthin <- terra::app(rlambda,fun=function(x)ifelse(is.na(x),eps,1-exp(-x)))
    out <- terra::mask(rthin,window)
  } else {
    out <- terra::mask(rlambda,window)
  }

  return(out)

}
