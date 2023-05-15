#' @title Non-parametric bootstrap to assess uncertainty in ppm parameters
#' @name bootstrap
#' @param object A ppm model object
#' @param nboot Default is 10; The number of bootstraps to run.
#' @param s If a specific value of lambda is required (maybe from the best model).
#' @param quiet Default is TRUE; Run the bootstrap with out any trace.
#' @param mc.cores Default is 1; the number of cores if run in parallel
#' @param \dots Ignored
#' @description This function uses a non-parametric bootstrap to estimate
#' uncertainty in a ppm. Conditional on X_j = {X_i,...,X_n}, let N*
#' have a Poisson distribution with Lambda_j. Lambda_j can be estimated as the
#' integral for each species intensities (lambda_ij), conditional on the
#' archetypes. For each species we can then draw \eqn{X*_{1},...,X*_{N}} by
#' sampling randomly with replacement N* times from X. We refit the ppm with
#' this new realisation of the ppm data. Estimated coefficients are returned
#' from each of B bootstrap model fits. This bootstrap object can then be used
#' in the to predict \link[stats]{predict} ppm and generate confidence intervals
#' function. If the model is very large (many po observations + quadrature sites)
#' vthis will be slow, and will take nboot times the time it takes to fit a
#' single ppm with the same data. So be warned for large datasets. We implement
#' method two from Cowling et al., (1996).
#' @references Cowling, A., Hall, P. and Phillips, M.J., 1996. Bootstrap
#' confidence regions for the intensity of a Poisson point process. Journal of
#' the American Statistical Association, 91(436), pp.1516-1524.
#' @importFrom stats rpois
#' @rdname bootstrap
#' @export bootstrap
#' @examples
#' \dontrun{
#'library(ppmData)
#'path <- system.file("extdata", package = "ppmData")
#'lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#'preds <- rast(lst)
#'presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
#'
#'ppmdata <- ppmData(npoints = 1000,presences=presences,
#'                   window = preds[[1]],
#'                   covariates = preds)
#'
#'form <- presence ~ poly(X,2) + poly(Y,2) + poly(max_temp_hottest_month,2)+
#'           poly(annual_mean_precip,2) + poly(annual_mean_temp,2) +
#'           poly(distance_from_main_roads,2)
#'
#'ft.ppm1 <- ppmFit(species_formula = form, ppmdata = ppmdata, method = 'lasso')
#'
#'bts <- bootstrap(ft.ppm1,nboot=100)
#'
#'## bootstraps for the 'best' model
#'cvl <- cvLambda(ft.ppm1)
#'
#'best_covars <- do.call(rbind,lapply(bts,function(x)x[,cvl$index[1]]))
#'
#'bts <- bootstrap(ft.ppm1, s = cvl$lambda.min)
#'}

"bootstrap" <- function (object, nboot=10, s=NULL, quiet=TRUE, mc.cores=1, ...){
  UseMethod("bootstrap", object)
}

#' @export
"bootstrap.ppmFit" <- function (object, nboot=10, s=NULL, quiet=TRUE, mc.cores=1, ...){

  if(is.null(object$titbits))
    stop("bootstrap needs `titbits=TRUE` when fitting a ppm.")

  fit0 <- predict(object = object)
  wtd_fitty <- fit0*object$ppmdata$ppmData$weights
  Lambda <- sum(wtd_fitty,na.rm=TRUE)

  coef_star <- plapply(X=1:nboot,
                       FUN=ppmBootFun,
                       Lambda=Lambda,
                       object=object,
                       s = s,
                       quiet=quiet,
                       .parallel = mc.cores,
                       .verbose = TRUE)

  # unlist(coef_star)
  boot.estis <- coef_star ## need a way to combine these at some point
  # class(boot.estis) <- "ppm_bootstrap"
  return(boot.estis)

}

"ppmBootFun" <- function(X, Lambda, object, s=NULL, quiet = TRUE){

  #simulate N from Lambda
  N_star <- rpois(n=1, lambda=Lambda)

  #find the scaled intensity for each of the observed spp locations
  dat_star <- object$ppmdata
  ids <- which(object$ppmdata$ppmData$presence==1)
  ids_star <- sample( ids, size=N_star, replace=TRUE)
  bootWts <- table( factor( ids_star, levels=ids))

  # if(object$titbits$fam=="binomial")
  #   bootWts <- ifelse(bootWts>0,1,0)

  #the bootstrap outcome
  dat_star$ppmData$presence[ids] <- bootWts

  #the bootstrap working variate (note that this extends the Berman and Turner as there are repeat presences in the bootstrap sample)
  # dat_star$ppmData$z <- dat_star$ppmData$presence / dat_star$ppmData$weights

  mf <- model.frame(formula = object$titbits$ppm_formula, data = dat_star$ppmData, weights = weights) # weights need to be name of weights in ppp
  mt <- terms(mf)
  X <- model.matrix(mt,mf)
  X <- X[,-1,drop=FALSE]
  y <- model.response(mf)

  wts <- model.weights(mf)

  offy <- model.offset(mf)
  if(is.null(offy))
    offy <- rep(0,length(y))

  if(is.null(s))
    s <- object$ppm$lambda


  if(object$titbits$method=="lasso"){
    ft <- glmnet::glmnet(x=X, y=y, weights = wts, offset = offy,
                         family = object$titbits$fam, alpha = 1,
                         lambda = s) #lasso
  }
  if(object$titbits$method=="ridge"){
    ft <- glmnet::glmnet(x=X, y=y, weights = wts, offset = offy,
                         family = object$titbits$fam, alpha = 0,
                         lambda = s) #ridge
  }

  boot.covars <- cbind(ft$a0,as.matrix(ft$beta))

  return(boot.covars)
}
