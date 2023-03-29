#'@title predict a ppmFit model
#'@name predict.ppmFit
#'@description This function should predict an intensity surface based on the
#'on the model fitting in the ppmFit.
#'@param object A fitted ppmFit object.
#'@param cvobject A cvLambda object. If NULL (default) the prediction function will run this internally.
#'@param newdata SpatRaster. A terra raster stack of covariates for the model, or it can be a data.frame.
#'@param type Character. Either "response","link" or "unit". The type of response variable to return. The default is 'response' which is on the intensity scale or 'link' which is one the linear predictor scale (log). Unit scales the intensity (response) by the area of each cell/prediction point.
#'@param offset Numeric vector or raster. If an offset is used in the model. Either an observed offset at prediction sites. If an offset is used and this is not known at prediction sites something like the mean offset used to fit the model can be used.
#'@param slambda Character Either 'lambda.min' or 'lambda.1se'. Value(s) of the penalty parameter lambda at which predictions are required. Default is "lambda.min".
#'@param quad.only Logical. If TRUE prediction is only done at the quadrature locations - useful for some of the diagnostic tools.
#'@param cores Integer. The number of cores to use in the prediction, useful for large rasters.
#'@param filename String Name of the raster file and path to save prediction. Default is NULL, otherwise it needs to be something like "pred.tif"
#'@param bigtif bool if TRUE it will try and do prediction via tiling, this will be slower but
#'will help with large tifs where holding the entire raster stack in memory is inpractical.
#'@param control list A list of control options for tiling. See the details below.
#'@param \\dots dots. Not used, but needed for prediction function.
#' @details For every large raster we can use tiling and parallel processing to do predictions
#' This is useful when making point process prediction for high resolution or broad scale rasters.
#' There a bunch of control arguments for running prediction using tiles can be passed as
#' follows to predict:
#' \describe{
#'  \item{predictionFile}{The path and file name of the prediction raster to be saved. Default is 'prediction.tif' and will be save in the current working directory.}
#'  \item{returnRaster}{Default is TRUE; leave this one alone unless you are using the predictWithTiles function on it's own}
#'  \item{tileFiles}{Name of the tile files, default is 'tile' and this will generate 'tile1.tif' to 'tileN.tif' where N is the total number of tiles to be generated.}
#'  \item{tilesDir}{The directory to store the tiles; default is a folder called 'tiles'}
#'  \item{vrtFile}{The vrt file needed to stitch together the files. Default is 'tmp.vrt'}
#'  \item{cacheTiles}{Should the covariate tiles created for prediction be cached? Default is FALSE. If TRUE the tiles used to generate the predictions will not be deleted once the function is finished. This means a user can reuse them for future predictions and skip the creation of tiles step.}
#'  \item{deleteTmp}{Will R delete all the tmp files created as part of the tiling? Default is TRUE. If false all input tiles, prediction tiles and vrt files will be retained. If cacheTiles is used then the input (covariate tiles) will be retained.}
#'  \item{ntiles}{The number of tiles to use across rows or columns, default is 10, which result in 100 tiles (10*10).}
#'  \item{mc.cores}{The number of cores to use when doing the tiles prediction. Default is 1, this will result in sequencial prediction per tile. If mc.cores > 1 tile prediction will be done in parallel.}
#'  }

#'@importFrom stats as.formula contrasts is.empty.model make.link predict rnorm runif sd var
#'@export
#'@examples
#'\dontrun{
#'library(ppmData)
#'library(ppmFit)
#'library(terra)
#'path <- system.file("extdata", package = "ppmData")
#'lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#'covariates <- rast(lst)
#'s <- sum(covariates)
#'covariates <- mask(covariates,s)
#'bias <- covariates[[1]]
#'names(bias) <- "bias"
#'covariates <- c(covariates,bias)
#'presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
#'
#'ppmdata <- ppmData(npoints = 10000,
#'                   presences=presences,
#'                   window = covariates[[1]],
#'                   covariates = covariates)
#'
#'sp_form <- presence ~ poly(annual_mean_precip,2,raw=TRUE) +
#'                      poly(annual_mean_temp,2,raw=TRUE) +
#'                      poly(distance_from_main_roads,2,raw=TRUE)
#'
#'## Fit a ppm using glmnet lasso
#' ft.ppm <- ppmFit(species_formula = sp_form, ppmdata=ppmdata)
#'
#' ## predict to the SpatRaster object
#' pred1 <- predict(ft.ppm, covariates)
#'
#' ## prediction using tiles (for large rasters)
#' pred2 <- predict(ft.ppm, covariates, bigtif=TRUE, control=list(mc.cores=2))
#'
#' ## predict to presence & quadrature sites
#' pred4 <- predict(ft.ppm)
#'
#' ## predict to just the quadrature sites
#' pred5 <- predict(ft.ppm, quad.only=TRUE)}

predict.ppmFit <- function(object,
                           cvobject = NULL,
                           newdata = NULL,
                           offset = NULL,
                           type = c("response","link"),
                           slambda= c("lambda.min","lambda.1se"),
                           # bias.correct = FALSE,
                           quad.only = TRUE,
                           filename = NULL,
                           bigtif = FALSE,
                           control = list(),
                           ...){

  ## set the controls - mainly for the tiles malarkey
  control <- setTilesControl(control)

  # newdata <- covariates
  type <- match.arg(type)
  model <- object$titbits$method# class(object[[1]])[1]
  slambda <- match.arg(slambda)
  object.mod <- object[[1]]


  ## if a glmnet cv object is missing then run this function
    if(is.null(cvobject)){
    cvfit <- glmnet::cv.glmnet(object$titbits$X,
                               object$titbits$y,
                               weights = as.numeric(object$titbits$wts),
                               alpha= ifelse(model=="lasso",1,0),
                               family = "poisson")
    }


  ## check if data is supplied
  if (is.null(newdata)) {
      newdata <- getPredQuad(object,quad.only)
      wts <- newdata$wts
      newdata <- newdata$X
  }

  ## check for offset
  if(is.null(offset)){
    offy <- getPredOffset(object = object,
                          newdata = newdata,
                          quad.only = quad.only)
  }


  if(any(isa(newdata,"SpatRaster"))){

      if(bigtif){
        pred <- predictWithTiles(newdata = newdata,
                                 model = cvfit,
                                 predfun = glmnetPredictFun,
                                 ntiles = control$ntiles,
                                 ppmfit = object,
                                 type=type,
                                 slambda = slambda,
                                 predictionFile = control$predictionFile,
                                 tileFiles = control$tileFiles,
                                 tilesDir = control$tilesDir,
                                 vrtFile = control$vrtFile,
                                 deleteTmp = control$deleteTmp,
                                 returnRaster = TRUE,
                                 mc.cores = control$mc.cores,
                                 cacheTiles =  control$cacheTiles,
                                 memfrac = control$memfrac,
                                 ...)
      } else {
        pred <- glmnetPredictFun(model = cvfit,
                                 newdata = newdata,
                                 ppmfit = object,
                                 type = type,
                                 offy = offy,
                                 slambda = slambda)
      }

    savePrediction(pred,filename)

  } else {
    ## Do prediction on a data.frame
    pred <- predict(object = cvfit, newx = newdata, type = type, s=slambda, newoffset = offy)
  }

  return(pred)

}


## Wrapper for predicting glmnet to terra rast
glmnetPredictFun <- function(model,
                             newdata,
                             ppmfit,
                             offy=NULL,
                             type = c("response","link"),
                             slambda = c("lambda.min","lambda.1se")) {

  slambda <- match.arg(slambda)
  type <- match.arg(type)

  if(missing(ppmfit))
    stop("ppmFit model is missing")
  if(missing(model))
    stop("glmnet cross validation object is missing")
  if(missing(newdata))
    stop("newdata is missing for predictions")

  if(any(isa(newdata,"SpatRaster"))){
    newdat2 <- terra::as.data.frame(newdata,xy=TRUE,na.rm=FALSE)
    xy <- newdat2[,1:2] ## coordinates for raster
    newdat2 <- newdat2[,-1:-2] ##data.frame without coordinates
    non.na.sites <- stats::complete.cases(newdat2)
    non.na.ids <- which(non.na.sites)
  }

  form2 <- ppmfit$titbits$ppm_formula
  form2[[2]] <- NULL
  new.mf <- stats::model.frame(form2,newdat2[non.na.ids,])
  mt <- stats::delete.response(ppmfit$titbits$terms)
  newx <- stats::model.matrix(mt,new.mf)
  offy <- stats::model.offset(new.mf)
  if(is.null(offy))
    offy <- rep(0,nrow(newx))

  preds <- predict(object = model, newx = newx, s = slambda, type = type, newoffset=offy)

  if(any(class(newdata)=="SpatRaster")){
    xy$preds <- newdat2[,1]
    xy$preds[non.na.ids] <- preds
    pred <- terra::rast(xy,type="xyz",crs=terra::crs(newdata))
  }
  return(pred)
}

# first pass at a predict with tiles approach.
predictWithTiles <-  function(newdata,
                              model,
                              predfun,
                              ntiles = 6,
                              predictionFile = "prediction.tif",
                              tileFiles = "tile",
                              tilesDir = "tiles",
                              vrtFile = "tmp.vrt",
                              deleteTmp = TRUE,
                              returnRaster = TRUE,
                              mc.cores = 1,
                              cacheTiles = FALSE,
                              memfrac = 0.8,
                              ...){

  ## create dir for making tiles
  if(!dir.exists(tilesDir))
    dir.create(tilesDir)

  ## Check to see if you want the prediction tiles to be cached
  ff <- paste0(tilesDir,"/",tileFiles,seq_len(ntiles*ntiles),".tif")
  if(!all(file.exists(ff))){
    x <- terra::rast(extent=terra::ext(newdata), ncols=ntiles, nrows=ntiles)
    ff <- terra::makeTiles(newdata, x, paste0(tilesDir,"/",tileFiles,".tif"),
                           overwrite=TRUE, memfrac=memfrac)
  }

  ## create the preds tmp files
  pff <- paste0(tilesDir,"/pred_",tileFiles,seq_len(ntiles*ntiles),".tif")

  predTile <- function(ii){
    allna <- all(is.na(terra::values(terra::rast(ff[ii]))))
    if(allna){
      r <- terra::rast(ff[ii])
      pred <- r[[1]]
      names(pred) <- "preds"
      terra::values(pred) <- NA
      terra::writeRaster(x = pred, filename = pff[ii], overwrite=TRUE)
    } else {
      tiledat <- terra::rast(ff[ii])
      pred <- predfun(model,
                      tiledat,
                      ...)
      names(pred) <- "preds"
      terra::writeRaster(x = pred, filename = pff[ii], overwrite=TRUE)
    }
  }

  # if(mc.cores > 1){
    message('Predicting tiles')
    plapply(seq_along(ff), function(ii){predTile(ii)}, .parallel = mc.cores, .verbose = TRUE)
  # } else {
    # for (ii in seq_along(ff)){
      # predTile(ii)
      # cat("predicted tile",ii,"of",ntiles^2,"\n")
    # }
  # }
  # }

  pred.merge <- terra::vrt(pff, paste0(tilesDir,"/",vrtFile), overwrite=TRUE)
  names(pred.merge) <- "prediction"

  if(!is.null(predictionFile)){
    terra::writeRaster(x = pred.merge,
                       filename = predictionFile,
                       overwrite=TRUE)
  }

  if(deleteTmp){

    unlink(pff)
    unlink(paste0(tilesDir,"/",vrtFile))

    if(!cacheTiles){
      unlink(ff)
      unlink(tilesDir)
    }

  }

  if(returnRaster){
    pred.out <- terra::rast(predictionFile)
    return(pred.out)
  }

}

# get the controls for tiles.
setTilesControl <- function(control){

  if (!("predictionFile" %in% names(control)))
    control$predictionFile <- "prediction.tif"
  if (!("returnRaster" %in% names(control)))
    control$returnRaster <- TRUE
  if (!("tileFiles" %in% names(control)))
    control$tileFiles <- "tile"
  if (!("tilesDir" %in% names(control)))
    control$tilesDir <- "tiles"
  if (!("vrtFile" %in% names(control)))
    control$vrtFile <- "tmp.vrt"
  if (!("cacheTiles" %in% names(control)))
    control$cacheTiles <- FALSE
  if (!("deleteTmp" %in% names(control)))
    control$deleteTmp <- TRUE
  if (!("ntiles" %in% names(control)))
    control$ntiles <- 10
  if (!("mc.cores" %in% names(control)))
    control$mc.cores <- 1
  if(!("memfrac" %in% names(control)))
    control$memfrac <- 0.8

  return(control)

}

getPredQuad <- function(object, quad.only){

  ## ppmlasso
  if(object$titbits$method=="ppmlasso"){
    if(quad.only){
      X <- as.data.frame(object$ppm$data[object$ppm$pres==0,])
      wts <- object$ppm$wt[object$ppm$pres==0]
    } else {
      X <- as.data.frame(object$ppm$data)
      wts <- object$ppm$wt
    }
  } else {
    if(quad.only){
      X <- object$titbits$X[object$titbits$y==0,]
      wts <- object$titbits$wts[object$titbits$y==0]
    } else {
      X <- object$titbits$X
      wts <- object$titbits$wts
    }
  }


  return(list(X=X,wts=wts))
}

getPredOffset <- function(object, newdata, quad.only){

                      if(!is.null(newdata)){
                        if(any(class(newdata)=="SpatRaster")){
                          offy <- newdata[[1]]
                          id <- !is.na(offy[])
                          roffy <- ifelse(id,0,NA)
                          offy <- terra::setValues(offy,roffy)
                          } else {
                          offy <- rep(0,nrow(newdata))
                        }
                      } else {
                        if(quad.only){
                            offy <- object$titbits$offy[object$titbits$y==0]
                          } else {
                            offy <- object$titbits$offy
                          }
                        }

  return(offy)
}

savePrediction <- function(pred,filename=NULL){
  if(any(class(pred)=="SpatRaster")){
    if(!is.null(filename))
      terra::writeRaster(x = pred, filename = filename, filetype ="GTiff", overwrite=TRUE)
  }
}

#'@title transform predictions from a ppmFit model.
#'@rdname transform
#'@name transform
#'@description Sometimes we want to transform the expectation/intensity of a point process
#'to something that scales between zero and one. Typically we might do this with a
#'logistic or complementary log-log transform. One challenge with this approach is that we often
#'don't know the true prevalence of a species in the landscape. We can thus transform the
#'intensity to a 'probability of presence' and can use the log(Lambda) as an constant offset based
#'the expected count. Ideally you need presence-absence data to working out this prevalence, but
#'this is typically missing.
#'@param object A fitted ppmFit model.
#'@param prediction A prediction from a ppmFit model, can either be a terra "SpatRaster" or data.frame/matrix
#'@param type What why of transform to do? The options are 'log', 'logit' and 'cloglog'. It assumes the input is the intensity prediction from the models.
#'

transform <- function(prediction, object, type= c("log","logit","cloglog")){

  type <- match.arg(type)

  if(isa(object,"SpatRaster")){
    if(type=="log"){
      pred.out <- log(prediction)
    }
    if(type=="logit"){
      Lambda <- terra::global(prediction*prod(terra::res(prediction)),'sum',na.rm=TRUE)
      pred.out <- 1/(1+exp(-log(prediction*prod(terra::res(prediction)))-log(as.numeric(Lambda))))
    }
    if (type == "cloglog"){
      Lambda <- terra::global(prediction*prod(terra::res(prediction)),'sum',na.rm=TRUE)
      pred.out <- 1-exp(-(exp(log(prediction*prod(terra::res(prediction)))+log(as.numeric(Lambda)))))
    }
  }else{
    stop("Currently only works as using a SpatRaster")
    # if(type=="log"){
    #   pred.out <- log(prediction)
    # }
    # if(type=="logit"){
    #   Lambda <- terra::global(prediction*prod(terra::res(prediction)),'sum',na.rm=TRUE)
    #   pred.out <- 1/(1+exp(-log(prediction*prod(terra::res(prediction)))-log(as.numeric(Lambda))))
    # }
    # if (type == "cloglog"){
    #   Lambda <- terra::global(prediction*prod(terra::res(prediction)),'sum',na.rm=TRUE)
    #   pred.out <- 1-exp(-(exp(log(prediction*prod(terra::res(prediction)))+log(as.numeric(Lambda)))))
    # }
  }

}


