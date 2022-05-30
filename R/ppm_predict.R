#'@title predict a ppmFit model
#'@name predict.ppmFit
#'@description This function should predict an intensity surface based on the
#'on the model fitting in the ppmFit.
#'@param object A fitted ppmFit object.
#'@param newdata SpatRaster. A terra raster stack of covariates for the model, or it can be a data.frame.
#'@param type Character. Either "response","link", "unit" & "cloglog". The type of response variable to return. The default is 'response' which is on the intensity scale or 'link' which is one the linear predictor scale (log).
#'@param offset Numeric vector or raster. If an offset is used in the model. Either an observed offset at prediction sites. If an offset is used and this is not known at prediction sites something like the mean offset used to fit the model can be used.
#'@param slambda Character Either 'lambda.min' or 'lambda.1se'. Value(s) of the penalty parameter lambda at which predictions are required. Default is "lambda.min".
#'@param quad.only Logical. If TRUE prediction is only done at the quadrature locations - useful for some of the diagnostic tools.
#'@param cores Integer. The number of cores to use in the prediction, useful for large rasters.
#'@param filename String Name of the raster file and path to save prediction. Default is NULL, otherwise it needs to be something like "pred.tif"
#'@param \\dots Additional parameter calls.
#'@export
#'@examples
#'library(ppmData)
#'library(ppmFit)
#'library(terra)
#'path <- system.file("extdata", package = "ppmData")
#'lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#'covariates <- rast(lst)
#'bias <- covariates[[1]]
#'names(bias) <- "bias"
#'covariates <- c(covariates,bias)
#'presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
#'ppmdata <- ppmData(npoints = 10000,presences=presences, window = covariates[[1]], covariates = covariates)
#'sp_form <- presence ~ poly(annual_mean_precip,2) + poly(annual_mean_temp,2) + poly(distance_from_main_roads,2)
#'## Fit a ppm using glmnet lasso
#' ft.ppm <- ppmFit(species_formula = sp_form, ppmdata=ppmdata)
#' pred <- predict(ft.ppm, covariates, type='cloglog')
#' pred <- predict(ft.ppm, type='response', quad.only=TRUE)


predict.ppmFit <- function(object,
                           # bootobject=NULL,
                           newdata = NULL,
                           type = c("response","link","unit","cloglog"),
                           offset = NULL,
                           # glmnet.cv = NULL,
                           slambda= c("lambda.min","lambda.1se"),
                           quad.only = TRUE,
                           cores = 1,
                           filename=NULL,
                           ...){

  # newdata <- covariates
  type <- match.arg(type)
  model <- object$titbits$method# class(object[[1]])[1]
  slambda <- match.arg(slambda)
  object.mod <- object[[1]]

  if(model%in%c("lasso","ridge")){
    # if(glmnet.cv){
    cvfit <- glmnet::cv.glmnet(object$titbits$x,
                               object$titbits$y,
                               weights = as.numeric(object$titbits$wts),
                               alpha= ifelse(model=="lasso",1,0),
                               family = "poisson")
    # }
  }

  ## check if data is supplied
  if (is.null(newdata)) {
      newdata <- getPredQuad(object.mod,quad.only)
      wts <- newdata$wts
      newdata <- newdata$X
      # wts <- newdata$wts
  }

  if(is.null(offset)){
    offset <- getPredOffset(object.mod = object,
                            newdata = newdata,
                            quad.only = quad.only)
  }


  if(any(class(newdata)=="SpatRaster")){
    ## if you pass a spatRaster stack
    ## let's predict as a raster
    ## Do predictions directly on the rasters with terra
    if(model=="ppmlasso")
      pred <- terra::predict(object = newdata, model=object.mod,
                             const = data.frame(wt = 1),
                             fun=pred.fun.ppmlasso, na.rm=TRUE, cores=cores)
    if(model=="glm")
      pred <- terra::predict(object = newdata, model=object.mod,
                             const = data.frame(weight = 1),
                             na.rm=TRUE, cores=cores)
    if(model=="gam")
      pred <- terra::predict(object = newdata, model=object.mod,
                             const = data.frame(weight = 1),
                             na.rm=TRUE, cores=cores)
    if(model=="lasso"){
      pred <- glmnetPredictFun(object = object,
                               cvobject = cvfit,
                               newdat = newdata,
                               offy = offy,
                               slambda = slambda)
    }
    if(model=="ridge"){
      pred <- glmnetPredictFun(object = object,
                               cvobject = cvfit,
                               newdat = newdata,
                               offy = offy,
                               slambda = slambda)
    }
    ## change the type
    if(type=="link")
      pred <- log(pred);
    if(type=="unit")
      pred <- pred*prod(terra::res(pred))
    if(type=="cloglog"){
      cell.pred <- pred*prod(terra::res(pred))
      Lambda <- terra::global(cell.pred,"sum",na.rm=TRUE)
      pred <- 1-exp(-pred/as.numeric(Lambda))
    }

    savePrediction(pred,filename)

  } else {
    ## Do prediction on a data.frame
    if(model=="ppmlasso")
      pred <- ppmlasso::predict.ppmlasso(object = object.mod, newdata = newdata)
    if(model=="glm")
      pred <- predict(object = object.mod, newdata = newdata, type = "response")
    if(model=="gam")
      pred <- predict(object = object.mod, newdata = newdata, type = "response")
    if(model=="lasso")
      pred <- predict(object = cvfit, newx = newdata, type = "response", s=slambda, newoffset = offy)
    if(model=="ridge")
      pred <- predict(object = cvfit, newx = newdata, type = "response", s=slambda, newoffset = offy)

    ## change the type
    if(type=="link")
      pred <- log(pred);
    if(type=="unit")
      pred <- pred*wts
    if(type=="cloglog"){
      tmp.pred <- pred*wts
      Lambda <- sum(tmp.pred)
      pred <- 1-exp(-pred/Lambda)
    }
  }
  pred
}

savePrediction <- function(pred,filename.raster=NULL){
  if(any(class(pred)=="SpatRaster")){
    if(!is.null(filename.raster))
      writeRaster(x = pred,filename = filename.raster,filetype ="GTiff",overwrite=TRUE)
  }
}

## Wrapper for predicting glmnet to terra rast
glmnetPredictFun <- function(object,
                             cvobject,
                             newdat,
                             offy=NULL,
                             slambda = c("lambda.min","lambda.1se")) {

  slambda <- match.arg(slambda)
  # type <- match.arg(type)

  if(missing(object))
    stop("ppmFit model is missing")
  if(missing(cvobject))
    stop("glmnet cross validation object is missing")
  if(missing(newdat))
    stop("newdat is missing for predictions")

  if(any(class(newdat)=="SpatRaster")){
    newdat2 <- terra::as.data.frame(newdat,xy=TRUE)
    xy <- newdat2[,1:2] ## cooridnates for raster
    newdat2 <- newdat2[,-1:-2] ##data.frame without coordinates
  }

  form2 <- object$titbits$ppm_formula
  form2[[2]] <- NULL
  new.mf <- stats::model.frame(form2,newdat2)
  mt <- stats::delete.response(object$titbits$terms)
  newx <- Matrix::sparse.model.matrix(mt,newdat2)

  if(is.null(offy))
    offy <- rep(0,nrow(newx))

  preds <- predict(object = cvobject, newx = newx, s = slambda, type = "response", newoffset=offy)

  if(any(class(newdat)=="SpatRaster")){
    pred <- rast(cbind(xy,preds),type="xyz")
  }

  return(pred)
}


ppmlassoPredictFun <- function(object, newdata=NULL, type = c("response","link"),
                              offset=NULL, interactions = NA, ...){

  if (any(lapply(newdata, class) == "factor")) {
    unpacknewdata = CatConvert(newdata)
    newdata = unpacknewdata$X
    cat.names = setdiff(unique(unpacknewdata$cat.names),
                        NA)
    use.form = as.character(object$formula)[2]
    for (i in 1:length(cat.names)) {
      use.form = gsub(cat.names[i], paste(names(newdata)[which(unpacknewdata$cat.names ==
                                                                 cat.names[i])], collapse = " + "), use.form)
    }
    object$formula = as.formula(paste("~", use.form))
  }
  var.0 = which(apply(newdata, 2, var) == 0)
  if (length(var.0) > 0) {
    newdata[, var.0] = newdata[, var.0] + rnorm(dim(newdata)[1],
                                                0, 1e-08)
  }
  mf = model.frame(object$formula, data = newdata)
  mt = attr(mf, "terms")
  X.des = if (!is.empty.model(mt))
    model.matrix(mt, mf, contrasts)
  else matrix(0, length(object$mu), 0L)
  X.var = X.des[, -1]
  if (is.null(object$s.means) == FALSE) {
    X.var = scale(X.var, center = object$s.means, scale = object$s.sds)
    X.des = cbind(1, X.var)
  }
  if (object$family == "area.inter") {
    if (is.na(interactions) == TRUE) {
      if (is.null(object$s.means) == FALSE) {
        X.des = cbind(X.des, min(scale(object$pt.interactions)))
      }
      if (is.null(object$s.means) == TRUE) {
        X.des = cbind(X.des, 0)
      }
    }
    if (is.na(interactions) == FALSE) {
      if (is.null(object$s.means) == FALSE) {
        X.des = cbind(X.des, scale(interactions, center = mean(object$pt.interactions),
                                   scale = sd(object$pt.interactions)))
      }
      if (is.null(object$s.means) == TRUE) {
        X.des = cbind(X.des, interactions)
      }
    }
  }

  link <- make.link('log')
  eta <- as.matrix(X.des) %*% object$beta + offy

  pred.int <- switch(type,
                     reponse = link$linkinv(eta),
                     link = eta)
  return(pred.int)

}

getPredQuad <- function(object, quad.only){

  ## ppmlasso
  if(any(class(object$ppm)=="ppmlasso")){
    if(quad.only){
      X <- as.data.frame(object$data[object.mod$pres==0,])
      wts <- object.mod$wt[object$pres==0]
    } else {
      X <- as.data.frame(object$data)
      wts <- object.mod$wt
    }
  }
  ## glmnet
  if(any(class(object$ppm)=="glmnet")){
    if(quad.only){
      X <- object$titbits$x[object$titbits$y==0,]
      wts <- object$titbits$wts[object$titbits$y==0]
    } else {
      X <- object$titbits$x
      wts <- object$titbits$wts
    }
  }


  return(list(X=X,wts=wts))
}

getPredOffset <- function(object.mod, newdata, quad.only){

                      if(!is.null(newdata)){
                        if(any(class(newdata)=="SpatRaster")){
                          offy <- newdata[[1]]
                          id <- !is.na(offy[])
                          roffy <- ifelse(id,0,NA)
                          offy <- setValues(offy,roffy)
                          } else {
                          offy <- rep(0,nrow(newdata))
                        }
                      } else {
                        if(quad.only){
                          if(class(object.mod)[1]=="ppmlasso"){
                            offy <- rep(0,nrow(object.mod$data[object.mod$pres==0,]))
                          }
                        } else {
                          if(class(object.mod)[1]=="ppmlasso"){
                            offy <- rep(0,nrow(object.mod$data))
                          }
                        }
                      }
  return(offy)
}





