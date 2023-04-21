# remotes::install_github("rspatial/terra")
library(terra)
library(sp)
df <- read.csv("dev/reptilia_2022-10-06.csv")
head(df)
ausProj <- "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
wgs84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
xy84 <- sp::SpatialPoints(df[,c("longitude","latitude")], proj4string=CRS(wgs84))
x <- spTransform(xy84, ausProj)
presences <- coordinates(x)
window <- rast("dev/test.tif")

effortOffset <- function(presences, window, fudge.factor=1e-9){

  v <- terra::vect(presences,crs=ausProj)

  rlambda <- terra::rasterize(v, window, fun=sum)

  rthin <- terra:::app(rlambda,fun=function(x)ifelse(is.na(x),fudge.factor,1-exp(-x)))

  rthin.out <- terra::mask(rthin,window)

  return(rthin.out)

}
