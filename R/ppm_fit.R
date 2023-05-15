#'@title A function for fitting various ppms using a ppmData data object
#'@name ppmFit
#'@rdname ppmFit
#'@description This is a wrapper function which will fit a ppm model using a range of existing fitting tools.
#'@param species_formula The ppm model formula for the environmental/niche relationships. What drives the distribution of the species?
#'@param bias_formula Default is NULL. The idea will be to implement this as part of an integrated model.
#'Currently its just a nice way to keep track of which covariates are used in species or bias variables.
#'If you include this formula at the moment it will be merged into the species formula and used additively in a loglinear equation
#'@param ppmdata A ppmData data object
#'@param family Character What family to use, current options are "poisson" for a poisson point process and "binomial" of an infinitely Weighted Logistic Regression (IWLR).
#'@param link Character What link function to use? Uses "default", "log" for poisson or "logit" for binomial. A user can specific "clogclog" for a binomial model.
#'@param method A method to fit the a ppm. Default is 'lasso', the alternative option is 'ridge' for ridge regression.
#'@param titbits A boolean call to assess if you want extra model objects to be returned.
#'@param standardise Do you want to standardise the data before modelling it?
#'@param control Options to pass to fitting functions.
#'@param \\dots Other things, not really used.
#'@author Skipton Woolley
#'@references Berman, M. and Turner, T.R., 1992. Approximating point process likelihoods with GLIM. Journal of the Royal Statistical Society: Series C (Applied Statistics), 41(1), pp.31-38.
#'@references Warton, D.I. and Shepherd, L.C., 2010. Poisson point process models solve the" pseudo-absence problem" for presence-only data in ecology. The Annals of Applied Statistics, pp.1383-1402. \url{https://doi.org/10.1214/10-AOAS331}
#'@references Renner, I.W. and Warton, D.I., 2013. Equivalence of MAXENT and Poisson point process models for species distribution modeling in ecology. Biometrics, 69(1), pp.274-281.
#'@details Uses the Berman-Turner device to fit an approximate loglike for PPM using a weighted Poisson model.
#'@importFrom stats binomial update model.frame terms model.matrix model.response model.offset model.weights weights update.formula poisson
#'@importFrom ppmData ppmData
#'@export
#'@examples
#'\dontrun{
#'library(ppmData)
#'library(terra)
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
#'ft.ppm1 <- ppmFit(species_formula = form, ppmdata=ppmdata, method='lasso')
#'ft.ppm2 <- ppmFit(species_formula = form, ppmdata=ppmdata, method='ridge')
#'ft.ppm3 <- ppmFit(species_formula = form, ppmdata=ppmdata, method='lasso',
#' family = 'binomial')
#'ft.ppm4 <- ppmFit(species_formula = form, ppmdata=ppmdata, method='ridge',
#' family = 'binomial', link="cloglog")
#'}

ppmFit <- function(species_formula = presence/weights ~ 1,
                   bias_formula = NULL,
                   ppmdata,
                   family = c("poisson","binomial"),
                   link =  c("default","log","logit","cloglog"),
                   method = c("lasso","ridge"),
                   control = list(),
                   titbits = TRUE,
                   standardise = FALSE,
                   ...){

  ## set control
  control <- setFitControl(control)

  # get the defaults
  family <- match.arg(family)
  link <- match.arg(link)
  method <- match.arg(method)

  if(standardise){

    Xscaled <- standardiseData(ppmdata)
    ppmdata$ppmData.scaled <- Xscaled

  }

  if(!isa(ppmdata,"ppmData"))
    stop("'ppmFit' requires a 'ppmData' object to run.")
  if(!any(c("lasso","ridge")%in%method)) #"ppmlasso"
    stop("'ppm.fit' must used one of the follow methods\n 'glm','gam','lasso',
         'ridge' to run.") #'ppmlasso'

  fam.out <- getFamily(family,link)
  fam <- fam.out$fam
  link <- fam.out$link

  ## merge formulas
  if(!is.null(bias_formula)){
   form <- merge.formula(species_formula,bias_formula)
  } else {
   form <- species_formula
  }

  ## check if response is ok and update depending on the family
  if(family=="poisson"){
    if(form[[2]]!="presence/weights"){
      form <- update(form,presence/weights ~ .)
    }
  }
  if(family=="binomial"){
    if(form[[2]]!="presence"){
      form <- update(form,presence ~ .)
    }
  }

  ## setup the model matrix/frames
  if(standardise){
    ppp <- ppmdata$ppmData.scaled
  } else {
    ppp <- ppmdata$ppmData
  }
  if(family=="binomial"){
    wts <- iwlrWeights(ppp)
    ppp$weights <- wts
  }

  mf <- model.frame(formula = form, data = ppp, weights = weights) # weights need to be name of weights in ppp
  mt <- terms(mf)
  X <- model.matrix(mt,mf)
  X <- X[,-1,drop=FALSE] ## drop intercept for glmnet
  y <- model.response(mf)

  ## get the weights - convert to IWLR if needed
  wts <- model.weights(mf)


  offy <- model.offset(mf)
  if(is.null(offy))
    offy <- rep(0,length(y))

  if(method=="lasso"){
    ft <- glmnet::glmnet(x=X, y=y, weights = wts, offset = offy,
                         family = fam, alpha = 1) #lasso
  }
  if(method=="ridge"){
    ft <- glmnet::glmnet(x=X, y=y, weights = wts, offset = offy,
                         family = fam, alpha = 0) #ridge
  }


  ## Save the titbits
  if(titbits){
    titbits <- list()
    titbits$species_formula <- species_formula
    titbits$bias_formula <- bias_formula
    titbits$ppm_formula <- form
    titbits$method <- method
    titbits$control <- control
    titbits$terms <- mt
    titbits$X <- X
    titbits$y <- y
    titbits$wts <- wts
    titbits$offy <- offy
    titbits$family <- family
    titbits$link <- link
    titbits$fam <- fam
  } else {
    titbits <- NULL
  }

  res <- list(ppm=ft,ppmdata=ppmdata,titbits=titbits)
  class(res)<- "ppmFit"
  return(res)

}

## check the families
getFamily <- function(family,link) {

    res <- list()
    if(family=="binomial" & link=="cloglog"){
      res$fam <- binomial(link="cloglog")
      res$link <- "cloglog"
    }
    else if(family=="binomial" & (link=="logit"|link=="default")){
      res$fam <- "binomial"
      res$link <- "logit"
    }
    else if(family=="poisson" & (link=="log"|link=="default")){
      res$fam <- "poisson"
      res$link <- "log"
    } else {
      stop (paste("'family'", family, "not recognised. Or combination of ",family," and ",link," not supported."))}
    return(res)

}

## check the weights
iwlrWeights <- function(ppmData){

  npres <- sum(ppmData$presence==1)
  nquad <- sum(ppmData$presence==0)
  wt <- ifelse(ppmData$presence == 1, 1, npres/nquad)

  return(wt)

}

# get the controls for tiles and other prediction things.
setFitControl <- function(control){

  if (!("mc.cores" %in% names(control)))
    control$mc.cores <- 1
  if (!("quiet" %in% names(control)))
    control$quiet <- FALSE

  return(control)

}



## function to standardize the covariate data
standardiseData <- function(ppmdata){

  Xscaled <- ppmdata$ppmData
  idx <- which(colnames(Xscaled)%in%c("ID","presence","weights"))
  Xscaled[,-idx,drop=FALSE] <- scale(Xscaled[,-idx,drop=FALSE])

  return(Xscaled)
}
