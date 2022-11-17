
#' @title cvfit to find the best penalty value (lambda) for glmnet ppm
#' @param object A ppm model object
#' @param \dots Other function calls for \link[glmnet]{cvfit}
#' @description This is just a wrapper for the \link[glmnet]{cvfit} function from the \link{glmnet} package.
#' @rdname cvLambda
#' @export cvLambda
#' @examples
#'library(ppmData)
#'library(terra)
#'library(glmnet)
#'path <- system.file("extdata", package = "ppmData")
#'lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#'preds <- rast(lst)
#'presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
#'
#'ppmdata <- ppmData(npoints = 1000,presences=presences,
#'                   window = preds[[1]],
#'                   covariates = preds)
#'
#'sp_form <- presence ~ poly(X,2) + poly(Y,2) + poly(max_temp_hottest_month,2)+
#'           poly(annual_mean_precip,2) + poly(annual_mean_temp,2) +
#'           poly(distance_from_main_roads,2)
#'
#'ft.ppm <- ppmFit(species_formula = sp_form, ppmdata=ppmdata, method='lasso')
#'
#'# Run the cross validation
#'cvl <- cvLambda(ft.ppm)
#'
#'# Plot the lambda cross-validation if so desired
#'plot(cvl)

cvLambda <- function (object, ...){
  UseMethod("cvLambda", object)
}

#' @export
cvLambda.ppmFit <- function(object, ...){

  if(is.null(object$titbits))
    stop("cvLambda needs titbits to run, please refit with 'titbits=TRUE'")

  cvfit <- glmnet::cv.glmnet(x = object$titbits$X,
                             y = object$titbits$y,
                             offset = object$titbits$offy,
                             weights = as.numeric(object$titbits$wts),
                             alpha= ifelse(object$titbits$method=="lasso",1,0),
                             family = object$titbits$fam,
                             ...)
  return(cvfit)

}

