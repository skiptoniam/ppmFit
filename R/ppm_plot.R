#' plotting function for a fitted weighted poisson
#' @name plot.ppmfit
#' @title Plotting a fitted spatial point process model.
#' @param model a fitted ppm model object.
#' @export
#' @examples

# library(ppmData)
# path <- system.file("extdata", package = "ppmData")
# lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
# preds <- terra::rast(lst)
# presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
#
# ## Plot study area/first raster and snails.
# plot(preds[[1]])
# points(presences,pch=16,cex=0.5,col=gray(0.4))
#
# ## Set up the raster stack (here I'm include coordinates as predictors)
# xy.st <- lonlat_from_window(preds[[1]],mask.na = TRUE)
# preds2 <- stack(xy.st,preds)
# ppmdata <- ppmData(npoints = 1000,presences=presences, window = preds[[1]], covariates = preds2)
#
# predict.ppmfit <- function(model,
#                         window,
#                         covarites,
#                         type = "response"){
#
#
#
# }


