#'@title A function for fitting various ppms using a ppmData data object
#'@name ppmFit
#'@rdname ppmFit
#'@description This is a wrapper function which will fit a ppm model using a range of existing fitting tools.
#'@param species_formula The ppm model formula for the environmental/niche relationships. What drives the distribution of the species?
#'@param bias_formula Default is NULL. The idea will be to implement this as part of an integrated model.
#'Currently its just a nice way to keep track of which covariates are used in species or bias variables.
#'If you include this formula at the moment it will be merged into the species formula and used additively in a loglinear equation
#'@param ppmdata A ppmData data object
#'@param method A method to fit the a ppm. Default is 'lasso'. Others options are: 'glm','gam','lasso','ridge' - removed ppmlasso for now
#'@param control Options to pass to fitting functions.
#'@param \\dots Other things, not really used.
#'@author Skipton Woolley
#'@references Berman, M. and Turner, T.R., 1992. Approximating point process likelihoods with GLIM. Journal of the Royal Statistical Society: Series C (Applied Statistics), 41(1), pp.31-38.
#'@references Warton, D.I. and Shepherd, L.C., 2010. Poisson point process models solve the" pseudo-absence problem" for presence-only data in ecology. The Annals of Applied Statistics, pp.1383-1402. \url{https://doi.org/10.1214/10-AOAS331}
#'@references Renner, I.W. and Warton, D.I., 2013. Equivalence of MAXENT and Poisson point process models for species distribution modeling in ecology. Biometrics, 69(1), pp.274-281.
#'#'@details Uses the Berman-Turner device to fit an approximate loglike for PPM using a weighted Poisson model.
#'@importFrom stats update model.frame terms model.matrix model.response model.offset update.formula poisson
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
#'sp_form <- presence ~ poly(X,2) + poly(Y,2) + poly(max_temp_hottest_month,2)+
#'           poly(annual_mean_precip,2) + poly(annual_mean_temp,2) +
#'           poly(distance_from_main_roads,2)
#'
#'ft.ppm1 <- ppmFit(species_formula = sp_form, ppmdata=ppmdata, method='glm')
#'ft.ppm2 <- ppmFit(species_formula = sp_form, ppmdata=ppmdata, method='gam')
#'ft.ppm3 <- ppmFit(species_formula = sp_form, ppmdata=ppmdata,
#' method='ppmlasso',control=list(n.fit=50))
#'ft.ppm4 <- ppmFit(species_formula = sp_form, ppmdata=ppmdata, method='lasso')
#'ft.ppm5 <- ppmFit(species_formula = sp_form, ppmdata=ppmdata, method='ridge')
#'}

ppmFit <- function(species_formula = presence/weights ~ 1,
                    bias_formula = NULL,
                    ppmdata,
                    method=c("lasso","glm","gam","ridge","ppmlasso"),
                    control=list(n.fit=20),...){

  # lambda will be a vector of
  method <- match.arg(method)

  if(!class(ppmdata)=="ppmData")
    stop("'ppmFit' requires a 'ppmData' object to run.")
  if(!any(c("glm","gam","lasso","ridge")%in%method)) #"ppmlasso"
    stop("'ppm.fit' must used one of the follow methods\n 'glm','gam','lasso',
         'ridge' to run.") #'ppmlasso'

  ## merge formulas
  if(!is.null(bias_formula)){
   form <- merge.formula(species_formula,bias_formula)
  } else {
   form <- species_formula
  }
  ## check if response is ok
  if(form[[2]]!="presence/weights"){
    form <- update(form,presence/weights ~ .)
  }

  ## setup the model matrix/frames
  if(any(!method%in%c("gam","ppmlasso"))){
    ## set up the data for ppm fit
    ppp <- ppmdata$ppmData
    mf <- model.frame(formula = form, data = ppp)
    mt <- terms(mf)
    x <- model.matrix(mt,mf)
    y <- model.response(mf)
    wts <- ppp$weights
    offy <- model.offset(mf)
    if(is.null(offy))
      offy <- rep(0,length(y))
    }

  if(method=="ppmlasso"){
    ## just need to do a little house keeping to make sure the data names and fomulas line up.
    form <- update.formula(form, NULL ~ .)
    dat <- ppmdata$ppmData
    colnames(dat)[which(colnames(dat)=="presence")] <- "Pres"
    colnames(dat)[which(colnames(dat)=="weights")] <- "wt"
  }

  # suppress warning because of the poisson is not int warning.
  if(method=="glm"){
    ft <- suppressWarnings(glm2::glm.fit2(x = x, y = y/wts, weights = wts,
                                          offset = offy, family = poisson()))
  }
  if(method=="gam"){
    ft <- suppressWarnings(mgcv::gam(formula = form, data = ppmdata$ppmData,
                                     weights = ppmdata$ppmData$weights,
                                     family = poisson()))
  }
  if(method=="ppmlasso"){
    ft <- suppressWarnings(ppmlasso::ppmlasso(formula = form, data = dat,
                                              n.fits = control$n.fit,
                                              family="poisson")) ## maybe could sub in ppmlasso
  }
  if(method=="lasso"){
    ft <- glmnet::glmnet(x=x, y=y/wts, weights = wts, offset = offy,
                         family = "poisson", alpha = 1) #lasso
  }
  if(method=="ridge"){
    ft <- glmnet::glmnet(x=x, y=y/wts, weights = wts, offset = offy,
                         family = "poisson", alpha = 0) #ridge
  }

  titbits <- list()
  titbits$species_formula <- species_formula
  titbits$bias_formula <- bias_formula
  titbits$ppm_formula <- form
  titbits$method <- method
  titbits$control <- control
  titbits$terms <- mt
  if(any(!method%in%c("gam","ppmlasso"))){
    titbits$x <- x
    titbits$y <- y
    titbits$wts <- wts
    titbits$offy <- offy
  }

  res <- list(ppm=ft,ppmdata=ppmdata,titbits=titbits)
  class(res)<- "ppmFit"
  return(res)

}


