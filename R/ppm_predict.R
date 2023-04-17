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
#'@param bias.values A named list with scalar per list object to make with specific values to make the bias covariates constant for prediction. For example, if distance_from_road is to be zero for prediction then the bias.values = list('distance_from_roads'=0), for multiple covariates it would look something like list('x1'=0,'x2'=10). If bias.values=NULL (by default) then no correction is made and the default values of the input data are used.
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
#' pred2 <- predict(ft.ppm, covariates, control=list(mc.cores=3))
#'
#' ## predict to presence & quadrature sites
#' pred3 <- predict(ft.ppm)
#'
#' ## predict to just the quadrature sites
#' pred4 <- predict(ft.ppm, quad.only=TRUE)}

predict.ppmFit <- function(object,
                           newdata = NULL,
                           cvobject = NULL,
                           offset = NULL,
                           type = c("response","link"),
                           slambda= c("lambda.min","lambda.1se"),
                           # bias.values = NULL,
                           quad.only = FALSE,
                           filename = NULL,
                           control = list(),
                           ...){

  ## set the controls - mainly for the tiles malarkey
  control <- setControl(control)

  type <- match.arg(type)
  method <- object$titbits$method# class(object[[1]])[1]
  slambda <- match.arg(slambda)
  object.mod <- object[[1]]


  ## if a glmnet cv object is missing then run this function
  if(is.null(cvobject)){
      cvobject <- cvLambda(object)
  }

  ## How do we do prediction with a preset dir of tiles?
  if(is.character(newdata)){
    newdata_tiles_path <- newdata
  }

  ## check if data is supplied
  if (is.null(newdata)) {
      newdata <- getPredQuad(object,quad.only)
      wts <- newdata$wts
      newdata <- newdata$X
  }

  ## check for offset
  if(is.null(offset)){
    if(!is.character(newdata)){ #only create offset if the data is not a tile set
      offy <- getPredOffset(object = object,
                            newdata = newdata,
                            quad.only = quad.only)
    } else {
      offy <- NULL
    }
  }

  if(any(isa(newdata,"SpatRaster"))){
    pred <- predictWithTerra(ppm = object,
                             cvppm = cvobject,
                             newdata = newdata,
                             type=type,
                             slambda = slambda,
                             control = control)

  } else if (is.character(newdata)) {
    pred <- predictWithTiles(ppm =
                             newdata_tiles_path = newdata,
                             model = cvobject,
                             predfun = glmnetPredictFun,
                             offset_tiles_path = offy,
                             type = type,
                             slambda = slambda,
                             control = control,
                             ...)

  } else {    ## Do prediction on a data.frame
    pred <- predict(object = cvobject, newx = newdata, type = type, s=slambda, newoffset = offy)
  }

  return(pred)

}

#'@name createPredTiles
#'@title Function to generate tiles for large predictions using terra
#'@param spatRasters A single or multiple layer rasters loaded using terra
#'@param ntiles The number of tiles to use for tiling in this will be squared so ntiles=10 equals 100 tiles in total.
#'@param tileNames The name to be given to each of the tiles e.g. "tile" would results in "tile1.tif" to "tileN.tif"
#'@param tilesDir The directory to save tiles
#'@export

createPredTiles <- function(spatRasters, ntiles = 10, tileNames="tile", tilesDir="tiles"){

  ## Check the tile directory exists
  if(!dir.exists(tilesDir))
    dir.create(tilesDir)

  ## Check to see if you want the prediction tiles to be cached
  ff <- paste0(tilesDir,"/",tileNames,seq_len(ntiles*ntiles),".tif")
  if(!all(file.exists(ff))){
    x <- terra::rast(extent=terra::ext(spatRasters), ncols=ntiles, nrows=ntiles)
    ff <- terra::makeTiles(spatRasters, x, paste0(tilesDir,"/",tileNames,".tif"),
                           overwrite=TRUE, memfrac=0.9, gdal=c("COMPRESS=LZW"))
  }
}

## Wrapper for predicting glmnet to terra rast
glmnetPredictFun <- function(ppmfit,
                             newdata,
                             cvppmfit = NULL,
                             offy = NULL,
                             type = c("response","link"),
                             slambda = c("lambda.min","lambda.1se")) {

  slambda <- match.arg(slambda)
  type <- match.arg(type)

  if(missing(ppmfit))
    stop("ppmFit object is missing")
  if(missing(newdata))
    stop("newdata is missing for predictions")

  if(any(isa(newdata,"SpatRaster"))){
    newdat2 <- terra::as.data.frame(newdata,xy=TRUE,na.rm=FALSE)
    xy <- newdat2[,1:2,drop=FALSE] ## coordinates for raster
    newdat2 <- newdat2[ , -1:-2, drop=FALSE] ##data.frame without coordinates
    non.na.sites <- stats::complete.cases(newdat2)
    non.na.ids <- which(non.na.sites)
  }

  ## drop NA sites
  newdat3 <- newdat2[non.na.ids,]

  ## do the dummy rows stuff
  addDummy <- checkPolyMat(newdat3)
  if(addDummy){
    dumdat <- addDummyRows(newdat3)
    newdat4 <- dumdat$X
    idx <- dumdat$idx
  } else {
    newdat4 <- newdat3
    idx <- NULL
  }

  ## do the model matrix stuff
  form2 <- ppmfit$titbits$ppm_formula
  form2[[2]] <- NULL
  new.mf <- stats::model.frame(form2,newdat4)
  mt <- stats::delete.response(ppmfit$titbits$terms)
  newx <- stats::model.matrix(mt,new.mf)

  if(addDummy)
    newx <- removeDummyRows(idx,newx)


  ## convert the offset to vector
  if(any(isa(offy,"SpatRaster"))){
    offy2 <- as.numeric(as.matrix(terra::as.data.frame(offy)))
      if(length(offy2)==nrow(newx)){
        offyin <- offy2
      } else {
        offyin <- rep(0,nrow(newx))
      }
  } else {
    if(is.null(offy)){
     offyin <- stats::model.offset(new.mf)
    }
    if(is.null(offyin)){
     offyin <- rep(0,nrow(newx))
    }
  }

  if(!is.null(cvppmfit)){
    preds <- predict(object = cvppmfit, newx = newx, s = slambda, type = type, newoffset=offyin)
  } else {
    preds.all <- predict(ppmfit, newx = newx, newoffset = offyin, type=type, exact=TRUE)
    preds <- apply(preds.all, 1, mean, na.rm=TRUE)
  }

  if(any(class(newdata)=="SpatRaster")){
    xy$preds <- newdat2[,1]
    xy$preds[non.na.ids] <- preds
    pred <- terra::rast(xy,type="xyz",crs=terra::crs(newdata))
  }

  return(pred)
}

checkPolyMat <- function(mat, degree=2){

  factors <- NULL
  for( i in 1:ncol(mat)){
      factors[i] <- is.factor(mat[,i])
  }

  rnk <- NULL
  for(j in seq_len(ncol(mat[,!factors]))){
    xbar <- mean(mat[,j])
    x <- mat[,j,drop=FALSE] - xbar
    X <- outer(as.numeric(as.matrix(x)), 0L:degree, "^")
    QR <- qr(X)
    rnk[j] <- QR$rank
  }

  out <- FALSE
  if(length(which(rnk<degree))>0)
    out <- TRUE

  return(out)

  }


addDummyRows <- function(x, degree = 2){

  factors <- NULL
  for( i in 1:ncol(x)){
    factors[i] <- is.factor(x[,i])
  }

  x[nrow(x) + 1:degree ,factors] <- unique(x[nrow(x),factors])
  x[nrow(x) + ((1:degree)-degree) ,!factors] <- rnorm(ncol(x[nrow(x),!factors])*degree,sd=1e-6)

  idx <- 1:degree

  return(list(X=x,idx=idx))

}

## remove the dummy columns
removeDummyRows <- function(idx, mm){

  idxs <- nrow(mm) + (idx-max(idx))
  mm.no.dummy <- mm[-idxs, ,drop=FALSE]

  dn1 <- attr(mm,"dimnames")[[1]][-idxs]
  dn2 <- attr(mm,"dimnames")[[2]]

  attr(mm.no.dummy,"dimnames") <-   list(dn1,dn2)
  attr(mm.no.dummy,"assign") <- attr(mm,"assign")
  attr(mm.no.dummy,"contrasts") <- attr(mm,"contrasts")

  return(mm.no.dummy)

}

glmnetPredictTerra <- function(model,
                               newdata,
                               cvmodel=NULL,
                               offy=NULL,
                               type = c("response","link"),
                               slambda = c("lambda.min","lambda.1se")) {

  slambda <- match.arg(slambda)
  type <- match.arg(type)

  form2 <- model$titbits$ppm_formula
  form2[[2]] <- NULL
  new.mf <- stats::model.frame(form2,newdata)
  mt <- stats::delete.response(model$titbits$terms)
  newx <- stats::model.matrix(mt,new.mf)
  offy <- stats::model.offset(new.mf)
  if(is.null(offy))
    offset <- rep(0,nrow(newx))

  if(!is.null(cvmodel)){
    result <- predict(cvmodel, newx = newx, s = slambda, newoffset = offset, type = type)
  } else {
    results.all <- predict(model$ppm, newx=newx, newoffset=offset, type=type, exact=TRUE)
    result <- apply(results.all,1,mean,na.rm=TRUE)
  }
  return(result)
}

predictWithTerra <- function(ppm,
                             cvppm=NULL,
                             newdata,
                             offy=NULL,
                             type = c("response","link"),
                             slambda = c("lambda.min","lambda.1se"),
                             control){


  pred <- terra::predict(object=newdata,
                         model=ppm,
                         fun=glmnetPredictTerra,
                         cvmodel=cvppm,
                         type=type,
                         slambda=slambda,
                         na.rm=TRUE,
                         cores=control$mc.cores,
                         filename=control$predictionFile,
                         overwrite=control$overwrite)

  return(pred)

}




#first pass at a predict with tiles approach.
predictWithTiles <-  function(ppm,
                              newdata_tiles_path,
                              cvppm = NULL,
                              offset_tiles_path = NULL,
                              type,
                              slambda,
                              control,
                              ...){

  ## create the preds tmp files
  ## covariate tiles
  tiles <- list.files(newdata_tiles_path,pattern="*.tif$")
  tiles <- reorderString(tiles,"tile")
  tile_paths <-  paste0(newdata_tiles_path,"/",tiles)

  ## offset tiles
  if(!is.null(offset_tiles_path)){
     offies <- list.files(offset_tiles_path)
     offies <- reorderString(offies,"tile")
     offy_paths <- paste0(offset_tiles_path,"/",offies)
  } else {
     offy_paths <- NULL
  }

  ## set up prediction tiles
  pff <- paste0(control$predsDir,"/pred_",tiles)
  predsDir <- control$predsDir

  ## create pred directory
  if(!dir.exists(predsDir))
    dir.create(predsDir)

  ## prediction function for tiles.
  predTile <- function(ii){
    allna <- all(is.na(terra::values(terra::rast(tile_paths[ii]))))
    if(allna){
      r <- terra::rast(tile_paths[ii])
      pred <- r[[1]]
      names(pred) <- "preds"
      terra::values(pred) <- NA
      terra::writeRaster(x = pred, filename = pff[ii], overwrite=TRUE)
    } else {
      tiledat <- terra::rast(tile_paths[ii])
      if(!is.null(offy_paths)){
        offy <- terra::rast(offy_paths[ii])
      } else{
        offy <- NULL
      }
      pred <- ppmFit:::glmnetPredictFun(ppmfit = ppm,
                                        newdata = tiledat,
                                        cvppmfit = cvppm,
                                        offy = offy,
                                        type = type,
                                        slambda = slambda)
      names(pred) <- "preds"
      terra::writeRaster(x = pred, filename = pff[ii], overwrite=TRUE)
    }
  }

  message('Predicting tiles')
  plapply(seq_along(tile_paths), function(ii){predTile(ii)}, .parallel = control$mc.cores, .verbose = TRUE)

  ## use the vrt function to stitch together the prediction rasters.
  pred.merge <- terra::vrt(pff, paste0(predsDir,"/",control$vrtFile), overwrite=TRUE)
  names(pred.merge) <- "prediction"

  ## Save the raster
  if(!is.null(control$filename)){
    terra::writeRaster(x = pred.merge,
                       filename = predictionFile,
                       overwrite=TRUE)
  }

  ## Delete the tmp files.
  if(control$deleteTmp){
    unlink(pff)
    unlink(paste0(predsDir,"/",vrtFile))
  }

  ## return the raster (default is true)
  if(control$returnRaster){
    pred.out <- terra::rast(predictionFile)
    return(pred.out)
  }

}

# get the controls for tiles.
setControl <- function(control){

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
  if (!("deleteTmp" %in% names(control)))
    control$deleteTmp <- TRUE
  if (!("predsDir" %in% names(control)))
    control$predsDir <- "preds"
  if (!("predictionFile" %in% names(control)))
    control$predictionFile <- "prediction.tif"
  if (!("mc.cores" %in% names(control)))
    control$mc.cores <- 1


  return(control)

}

getPredQuad <- function(object, quad.only){

    if(quad.only){
      X <- object$titbits$X[object$titbits$y==0,,drop=FALSE]
      wts <- object$titbits$wts[object$titbits$y==0]
    } else {
      X <- object$titbits$X
      wts <- object$titbits$wts
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

reorderString <- function(x, split){

  split1 <- strsplit(x, "tile")
  split2 <- as.numeric(sapply(split1, function(y) y <- sub(".tif", "", y[2])))
  reorderedNames <- x[order(split2)]
  return(reorderedNames)

}

# savePrediction <- function(pred,filename=NULL){
#   if(any(class(pred)=="SpatRaster")){
#     if(!is.null(filename))
#       terra::writeRaster(x = pred, filename = filename, filetype ="GTiff", overwrite=TRUE)
#   }
# }

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
#'@export

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


