#'@title predict a ppmFit model
#'@name predict.ppmFit
#'@description This function should predict an intensity surface based on the
#'on the model fitting in the ppmFit.
#'@param object A fitted ppmFit object as fitted using \link{ppmFit}[ppmFit].
#'@param cvobject A cvLambda object. If NULL (default) the prediction will be calculated as the average of all lambda predictions.
#'@param newdata Either a SpatRaster stack with the same data as those in the model; a character string which points to a path which contains tiles as created using \link{ppmFit}[createPredTiles] or a data.frame
#'@param type Character. Either "response","link" or "unit". The type of response variable to return. The default is 'response' which is on the intensity scale or 'link' which is one the linear predictor scale (log). Unit scales the intensity (response) by the area of each cell/prediction point.
#'@param offset Either a SpatRaster with an offset that will be used in the model; a character string which points to a path which contains a tiled version of the offset created using \link{ppmFit}[createPredTiles] or numeric vector.
#'@param slambda Character Either 'lambda.min' or 'lambda.1se'. Value(s) of the penalty parameter lambda at which predictions are required. Default is "lambda.min".
#'@param quad.only Logical. If TRUE prediction is only done at the quadrature locations - useful for some of the diagnostic tools. Only works if cvobject is passed to the prediction function.
#'@param filename String A filename or complete file path and name when saving a raster prediction
#'@param control list A list of control options for tiling. See the details below.
#'@param \\dots dots. Not used, but needed for prediction function.
#' @details For very large rasters we can use tiling and parallel processing to do predictions.
#' This is for making prediction for high resolution or at broad spatial scales.
#' There a bunch of control arguments for running prediction using tiles can be passed as
#' follows to predict:
#' \describe{
#'  \item{predictionFile}{The path and file name of the prediction raster to be saved. Default is 'prediction.tif' and will be save in the current working directory.}
#'  \item{returnRaster}{Default is TRUE; leave this one alone unless you are using the predictWithTiles function on it's own}
#'  \item{tileFiles}{Name of the tile files, default is 'tile' and this will generate 'tile1.tif' to 'tileN.tif' where N is the total number of tiles to be generated.}
#'  \item{tilesDir}{The directory to store the tiles; default is a folder called 'tiles'}
#'  \item{vrtFile}{The vrt file needed to stitch together the files. Default is 'tmp.vrt'}
#'  \item{deleteTmp}{Should all the tmp files created as part of the tiling? Default is TRUE, can turn this off for trouble shooting}
#'  \item{predsDir}{The directory to hold the tmp prediction tiles. Default is "preds"}
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
  control <- setPredControl(control)

  ## Check for the filename and if not null overwrite the prediction path
  if(!is.null(filename)){
    control$predictionFile <- filename
  }

  type <- match.arg(type)
  method <- object$titbits$method# class(object[[1]])[1]
  slambda <- match.arg(slambda)
  object.mod <- object[[1]]

  ## check if data is supplied
  if (is.null(newdata)) {
      newdata <- getPredQuad(object,quad.only)
      wts <- newdata$wts
      newdata <- newdata$X
  }

  ## How do we do prediction with a preset dir of tiles?
  if(is.character(newdata)){
    newdata_tiles_path <- newdata
  }

  ## check for offset
  offy <- getPredOffset(offset = offset,
                        object = object,
                        newdata = newdata,
                        quad.only = quad.only)

  # print(offy)

  if(any(isa(newdata,"SpatRaster"))){
    pred <- predictWithTerra(ppm = object,
                             cvppm = cvobject,
                             newdata = newdata,
                             type=type,
                             slambda = slambda,
                             control = control)

  } else if (is.character(newdata)) {
    pred <- predictWithTiles(ppm = object,
                             cvppm = cvobject,
                             newdata_tiles_path = newdata,
                             offset_tiles_path = offy,
                             type = type,
                             slambda = slambda,
                             control = control)

  } else {    ## Do prediction on a data.frame

    if(!is.null(cvobject)){
      pred <- predict(object = cvobject, newx = newdata, s = slambda, type = type, newoffset = offy)
    } else {
      preds.all <- predict(object = object$ppm, newx = newdata, type=type,  newoffset = offy, exact=TRUE)
      pred <- apply(preds.all, 1, mean, na.rm=TRUE)
    }

  }

  return(pred)

}

#'@name createPredTiles
#'@title Function to generate tiles for large predictions using terra
#'@param spatRasters A single or multiple layer rasters loaded using terra
#'@param ntiles The number of tiles to use for tiling in this will be squared so ntiles=10 equals 100 tiles in total.
#'@param tileNames The name to be given to each of the tiles e.g. "tile" would results in "tile1.tif" to "tileN.tif"
#'@param tilesDir The directory to save tiles
#'@param overwrite Default is FALSE, should you delete and overwrite existing tiles?
#'@export

createPredTiles <- function(spatRasters, ntiles = 10, tileNames="tile", tilesDir="tiles", overwrite=FALSE){

  ## Check the tile directory exists
  if(!dir.exists(tilesDir))
    dir.create(tilesDir)

  ## Check to see if you want the prediction tiles to be cached
  ff <- paste0(tilesDir,"/",tileNames,seq_len(ntiles*ntiles),".tif")
  if(overwrite){
    unlink(ff)
  }
  if(!all(file.exists(ff))){
    x <- terra::rast(extent=terra::ext(spatRasters), ncols=ntiles, nrows=ntiles)
    ff <- terra::makeTiles(spatRasters, x, paste0(tilesDir,"/",tileNames,".tif"),
                           overwrite=overwrite, memfrac=0.9, gdal=c("COMPRESS=LZW"))
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
    newdata2 <- newdat2[ , -1:-2, drop=FALSE] ##data.frame without coordinates
    non.na.sites <- stats::complete.cases(newdat2)
    non.na.ids <- which(non.na.sites)

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
  } else {
    newdat4 <- newdata
  }
  # print(nrow(newdat4))


  ## do the model matrix stuff
  form2 <- ppmfit$titbits$ppm_formula
  form2[[2]] <- NULL
  new.mf <- stats::model.frame(form2,newdat4)
  mt <- stats::delete.response(ppmfit$titbits$terms)
  newx <- stats::model.matrix(mt,new.mf)

  if(any(isa(newdata,"SpatRaster"))){
    if(addDummy)
      newx <- removeDummyRows(idx,newx)
  }

  if(is.null(offy)){
    offy <- rep(0,nrow(newx))
  }


  ## convert the offset to vector
  if(any(isa(offy,"SpatRaster"))){
    offy2 <- terra::as.data.frame(offy,xy=TRUE,na.rm=FALSE)
    xyoffy <- offy2[,1:2,drop=FALSE] ## coordinates for raster
    offy2 <- offy2[ , -1:-2,drop=FALSE] ##data.frame without coordinates
    non.na.sites.offy <- stats::complete.cases(offy2)
    non.na.ids.offy <- which(non.na.sites)
    offy2 <- offy2[non.na.ids.offy,]

      if(length(offy2)==nrow(newx)){
        offy <- offy2
      } else {
        offy <- rep(0,nrow(newx))
      }
  }

  if(!is.null(cvppmfit)){
    preds <- predict(object = cvppmfit, newx = newx, s = slambda, type = type, newoffset = offy)
  } else {
    preds.all <- predict(object = ppmfit$ppm, newx = newx, type=type,  newoffset = offy, exact=TRUE)
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
  if(length(which(rnk<=degree))>0)
    out <- TRUE

  return(out)

  }


addDummyRows <- function(x, degree = 2){

  factors <- NULL
  for( i in 1:ncol(x)){
    factors[i] <- is.factor(x[,i])
  }

  x[nrow(x) + 1:degree ,factors] <- unique(x[nrow(x),factors])
  x[nrow(x) + ((1:degree)-degree) ,!factors] <- rnorm(n = ncol(x[,!factors,drop=FALSE])*degree,
                                                      mean = rep(apply(x[,!factors,drop=FALSE],2,mean,na.rm=TRUE),each=degree),
                                                      sd =rep(apply(x[,!factors,drop=FALSE],2,sd,na.rm=TRUE),each=degree)+1e-6)

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


predictWithTiles <-  function(ppm,
                              cvppm = NULL,
                              newdata_tiles_path,
                              offset_tiles_path = NULL,
                              type,
                              slambda,
                              control){

  ## create the preds tmp files
  ## covariate tiles
  tiles <- list.files(newdata_tiles_path,pattern="*.tif$")
  tiles <- reorderString(tiles,"tile")
  tile_paths <-  paste0(newdata_tiles_path,"/",tiles)

  ## offset tiles
  if(!is.null(offset_tiles_path)){
     offies <- list.files(offset_tiles_path,pattern="*.tif$")
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
      pred <- glmnetPredictFun(ppmfit = ppm,
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
  if(!is.null(control$predictionFile)){
    terra::writeRaster(x = pred.merge,
                       filename = control$predictionFile,
                       overwrite=TRUE)
  }

  ## Delete the tmp files.
  if(control$deleteTmp){
    unlink(pff)
    unlink(paste0(predsDir,"/",control$vrtFile))
  }

  ## return the raster (default is true)
  if(control$returnRaster){
    pred.out <- terra::rast(control$predictionFile)
    return(pred.out)
  }

}

# get the controls for tiles and other prediction things.
setPredControl <- function(control){

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
  if (!("overwrite" %in% names(control)))
    control$overwrite <- TRUE

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

getPredOffset <- function(offset, object, newdata, quad.only){

                  if(is.null(offset)){
                    if(!is.null(newdata)){
                      if(is.character(newdata)){
                        offy <- NULL
                      } else if(any(class(newdata)=="SpatRaster")){
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
                  } else {
                    offy <- offset
                  }


  return(offy)
}

reorderString <- function(x, split){

  split1 <- strsplit(x, "tile")
  split2 <- as.numeric(sapply(split1, function(y) y <- sub(".tif", "", y[2])))
  reorderedNames <- x[order(split2)]
  return(reorderedNames)

}

#'@title transform predictions from a ppmFit model.
#'@rdname transformPrediction
#'@name transformPrediction
#'@description Sometimes we want to transform the expectation/intensity of a point process
#'to something that scales between zero and one. Typically we might do this with a
#'logistic or complementary log-log transform. One challenge with this approach is that we often
#'don't know the true prevalence of a species in the landscape. We can thus transform the
#'intensity to a 'probability of presence' and can use the log(Lambda) as an constant offset based
#'the expected count. Ideally you need presence-absence data to working out this prevalence, but
#'this is typically missing.
#'@param prediction A prediction from a ppmFit model, can either be a terra "SpatRaster" or data.frame/matrix
#'@param type What why of transform to do? The options are 'log', 'logit' and 'cloglog'. It assumes the input is the intensity prediction from the models.
#'@export

transformPrediction <- function(prediction, type= c("log","logit","cloglog")){

  type <- match.arg(type)

  if(isa(prediction,"SpatRaster")){
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

  return(pred.out)

}


