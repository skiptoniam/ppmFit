#'@title A function for doing a very basic cross validation fitting based on a ppmData data object
#'@name cvFit
#'@rdname cvFit
#'@description This is a wrapper function which will fit a ppm model using a range of existing ppm fitting tools.
#'@param species_formula The ppm model formula.
#'@param bias_formula Default is NULL. The idea will be to implement this as part of an integrated model.
#'Currently its just a nice way to keep track of which covariates are used in species or bias variables.
#'If you include this formula at the moment it will be merged into the species formula.
#'@param ppmdata ppmData object
#'@param family Character What family to use, current options are "poisson" for a poisson point process and "binomial" of an infinitely Weighted Logistic Regression (IWLR).
#'@param link Character What link function to use? Uses "default", "log" for poisson or "logit" for binomial. A user can specific "clogclog" for a binomial model.
#'@param method A method to fit the a ppm. Default is 'lasso', the alternative option is 'ridge' for ridge regression.
#'@param titbits A boolean call to assess if you want extra model objects to be returned.
#'@param type character Use a "thin" approch or a "block" approach for setting up a cross-validation
#'@param control list to pass to fitting functions.
#'@param nsim integer The number of cross validations to do (default is five).
#'@param p numeric The thinning probability (0,1).
#'@param res numeric The resolution of blocks if you want to try block resample
#'@param seed integer An integer to uses for set.seed(seed)
#'@param \\dots other things, not really used.
#'@importFrom utils txtProgressBar setTxtProgressBar
#'@author Skipton Woolley
#'@export
#'@examples
#'\dontrun{
#'library(ppmData)
#'library(terra)
#'path <- system.file("extdata", package = "ppmData")
#'lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
#'preds <- rast(lst)
#'presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
#'ppmdat <- ppmData(npoints = 10000,presences=presences, window = preds[[1]],
#' covariates = preds)
#'sp_form <- presence ~ poly(X,2) + poly(Y,2) +
#' poly(max_temp_hottest_month,2) + poly(annual_mean_precip,2) +
#' poly(annual_mean_temp,2) + poly(distance_from_main_roads,2)
#'cvft.ppm <- cvFit(species_formula = sp_form, ppmdata=ppmdat, family = "binomial", link = "logit")
#'}
cvFit <- function(species_formula = presence/weights ~ 1,
                  bias_formula = NULL,
                  ppmdata,
                  family = c("poisson","binomial"),
                  link =  c("default","log","logit","cloglog"),
                  method = c("lasso","ridge"),
                  control = list(),
                  titbits = TRUE,
                  type=c("thin","block"),
                  nsim=5,
                  p=0.25,
                  res=1,
                  seed=NULL,
                  ...){

    ## set control
    control <- setFitControl(control)

    ## set up a default method
    method <- match.arg(method)
    type <- match.arg(type)
    family <- match.arg(family)
    link <- match.arg(link)

    fam.out <- getFamily(family,link)
    fam <- fam.out$fam
    link <- fam.out$link

    #need update the quadrature per-fit
    cvdatasets <- list()
    for(ii in seq_len(nsim)){
      if(type=="thin")
        cvdatasets[[ii]] <- pThin(ppmdata,p=p)
      if(type=="block")
        cvdatasets[[ii]] <- blockSample(ppmdata, p=p, res=res, seed = seed)
    }


    if( !control$quiet&&nsim>1){
      # if(nsim>1){
      chars <- c("_@_'' ","><((('>")
      pb <- utils::txtProgressBar(min = 1, max = nsim, style = 3, char = chars[sample(2,1)])
      # }
    }

    fit.fun <- function(ii){
      if( !control$quiet&&nsim>1)
        utils::setTxtProgressBar(pb, ii)
        ppmfit <- ppmFit(species_formula = species_formula,
                            bias_formula = bias_formula,
                            ppmdata = cvdatasets[[ii]]$train,
                            family = family,
                            method = method,
                            control = control)
      return(ppmfit)
    }

    cvppmfits <- list()
    for(ii in seq_len(nsim)){
      cvppmfits[[ii]] <- fit.fun(ii)
    }

    res <- list(cvfits=cvppmfits,cvdata=cvdatasets)
    class(res) <- "cvFit"
    return(res)
}

#' @title A function for p-thinning a point process
#' @name pThin
#' @rdname pThin
#' @description We need to run a monte-carlo/kfold cross validation on the observed
#' point process. To do so requires the thinning of the observed pp. We can do
#' this via independent p-thinning. If the thinning is independent then the
#' resulting thinned point process can be used as another point process!
#' The current realization we have is the observed species locations. So if we
#' thin this based on a p-thinning and use the retained points as the validation
#' set and the thinned (removed) as the testing set we should be right.
#' @param object ppmData object that represents our species point process
#' @param p double The thinning probability (Bernoulli variable).
#' @param seed integer A seed if you want to set seed for randomisation
#' @author Skipton Woolley

pThin <- function(object, p = 0.25, seed=NULL){

   ## get presence IDs
   z <- object$ppmData$presence
   sites <- object$ppmData[z==1,object$params$coord]

   ## thin sites
   cvsites <- pthin(sites, p = p, seed = seed, ids = TRUE)

   ## create two new data objects
   # z_pres <- z[z==1]
   z_quad <- which(z==0)
   z_train <- which(cvsites$retained_sites)
   z_test <-  which(cvsites$thinned_sites)

   ppmdat_train <- ppmdat_test <- object
   # ppmdat_test <- ppmdat

   ppmdat_train$ppmData <- object$ppmData[c(z_train,z_quad),]
   ppmdat_test$ppmData <- object$ppmData[c(z_test,z_quad),]

   # ppmdat_train$presences

   res <- list()
   res$train <- ppmdat_train
   res$test <- ppmdat_test
   return(res)

}

pthin <- function(sites, p=0.25, seed = NULL, ids = FALSE){


  if(!ncol(sites)==2)stop("'sites' should just be the coordinates of the species presences.")

  # The number of
  n <- nrow(sites)

  # Bernoulli variable of known points being thinned
  if(!is.null(seed))set.seed(seed)

  boole_thin <- runif(n)<p;
  boole_retain <- !boole_thin

  if(ids){
    res <- list(thinned_sites=boole_thin,
                retained_sites=boole_retain)
    return(res)
  }


  #locations of thinned points
  thinned=sites[boole_thin,]
  retained=sites[boole_retain,]

  # return as list
  res <- list(thinned_sites=thinned,
              retained_sites=retained)

  return(res)

}

#' @title block sample of ppmData object
#' @name blockSample
#' @description Do a block sample of sites if you want a tougher k-folds test
#' @param object ppmData A ppmData object
#' @param p numeric The thinning probability (Bernoulli variable); values are (0,1).
#' @param res numeric The resolution (size in map coordinates) of the blocks.
#' @param seed numeric The seed to set if you want to set.seed(). Default is NULL (no seed)
#' @author Skipton Woolley

blockSample <- function(object, p = 0.25, res = 1, seed=NULL){

  ## get presence IDs
  # z <- which(object$ppmData$presence==1)
  sites <- object$ppmData[,object$params$coord]

  ## block thin sites
  block_ids <- blocksample(sites = sites,
                         coords = object$params$coord,
                         p.blocks = p,
                         block.res = res,
                         seed = seed,
                         ids = TRUE)

  # plot(cvsites$block_train_sites[!rownames(cvsites$block_train_sites)%in%z,],pch=19,cex=0.5)
  # points(cvsites$block_train_sites[rownames(cvsites$block_train_sites)%in%z,],pch=19,cex=0.5,col='red')

  # plot(cvsites$block_test_sites[!rownames(cvsites$block_train_sites)%in%z,],pch=19,cex=0.5)
  # points(cvsites$block_test_sites[rownames(cvsites$block_train_sites)%in%z,],pch=19,cex=0.5,col='red')

  ppmdat_train <- ppmdat_test <- object
  ppmdat_train$ppmData <- object$ppmData[c(block_ids$block_train),]
  ppmdat_test$ppmData <- object$ppmData[c(block_ids$block_test),]

  res <- list()
  res$train <- ppmdat_train
  res$test <- ppmdat_test
  return(res)

}


blocksample <- function(sites, coords=c("X","Y"), p.blocks = 0.25, block.res = 1, seed = NULL, ids=FALSE){


  cell.group <- rep(0, length(sites[,coords[1]]))
  n.groups   <- ceiling(c(max((sites[,coords[1]] - min(sites[,coords[1]]))/block.res), max((sites[,coords[2]] - min(sites[,coords[2]]))/block.res)))

  xq <- block.res*(1:n.groups[1]) + min(sites[,coords[1]])
  for (i.group in 1:n.groups[1]){
    cell.group <- cell.group + as.numeric(sites[,coords[1]] > xq[i.group])
  }

  yq <- block.res*(1:n.groups[2]) + min(sites[,coords[2]])
  for (i.group in 1:n.groups[2]){
    cell.group <- cell.group + n.groups[1] * as.numeric(sites[,coords[2]] > yq[i.group])
  }

  block.group <- factor(cell.group)

  cat("There are a total of",length(levels(block.group)), "blocks at a",block.res,"resolution")

  if(!is.null(seed))set.seed(seed)

  levs  <- levels(block.group)
  nlevs <- length(levs)
  p.blocks.in <- ceiling(p.blocks*nlevs)
  lev.sample <- sample(levs,p.blocks.in)

  block.sites <- which(block.group%in%lev.sample)

  if(ids){
    res <- list(block_train=which(!block.group%in%lev.sample),
                block_test=which(block.group%in%lev.sample))
    return(res)
  }

  #locations of thinned points
  block_train=sites[-block.sites,]
  block_test=sites[block.sites,]

  # return as list
  res <- list(block_train_sites=block_train,
              block_test_sites=block_test)

  return(res)
}

## function to do AUC/ROC ect.
#' @title Assess performance of models using hold-out dataset
#' @param object A cvFit model object, which is essentially a bunch of model fitted to a holdout dataset.
#' @param slambda A numeric value to represent 'lambda.min', 'lambda.1se' or a set of lambda valies. Value(s) of the penalty parameter lambda at which predictions are required. Default is "lambda.min".
#' @param \dots Ignored
#' @description This function produces summary performance measures for the glmnet model(s). It requires a test dataset.
#' @name cvTests
#' @rdname cvTests
#' @export cvTests
"cvTests" <- function (object, slambda, ...){
  UseMethod("tests", object)
}

#' @export
"cvTests" <- function(object, slambda= NULL, ...){

  if(!isa(object,"cvFit"))
    stop("cvTests needs a cvFit object to work.")

  ## get the number of cv fits
  ncv <- length(object$cvfits)

  cv.assess <- list()

  ## need a cvobject per fit
  # if(!is.null(cvobject)){
    # idmin = match($lambda.min, fit2c$lambda)
  # }

  for(ii in seq_len(ncv)){

    mf <- model.frame(formula = object$cvfits[[ii]]$titbits$ppm_formula,
                      data = object$cvdata[[ii]]$test$ppmData,
                      weights = weights) # weights need to be name of weights in ppp
    mt <- terms(mf)
    xnew <- model.matrix(mt,mf)
    xnew <- delete.intercept(xnew)
    ynew <- model.response(mf)
    wts <- model.weights(mf)
    offy <- model.offset(mf)
    if(is.null(offy))
      offy <- rep(0,length(ynew))

    cv.assess[[ii]] <- glmnet::assess.glmnet(object = object$cvfits[[ii]]$ppm,
                                           newx = xnew,
                                           newy = ynew,
                                           family = object$cvfits[[ii]]$titbits$family,
                                           weights = wts,
                                           newoffset = offy,
                                           s = slambda)


  }

  return(cv.assess)

}



