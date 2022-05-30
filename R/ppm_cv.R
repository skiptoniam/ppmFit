#'@title A function for doing a very basic cross validation fitting based on a ppmData data object
#'@name cvfit
#'@rdname cvfit
#'@description This is a wrapper function which will fit a ppm model using a range of existing ppm fitting tools.
#'@param species_formula The ppm model formula.
#'@param bias_formula Default is NULL. The idea will be to implement this as part of an integrated model.
#'Currently its just a nice way to keep track of which covariates are used in species or bias variables.
#'If you include this formula at the moment it will be merged into the species formula.
#'@param object A ppmData object
#'@param method A method to fit the a ppm. Default is 'glm'. Others options are: 'gam','ppmlasso','lasso','ridge'
#'@param control Options to pass to fitting functions.
#'@param nsim Integer The number of cross validations to do (default is five).
#'@param ncores Integer The number of cores to used for parallel computing
#'@param \\dots Other things, not really used.
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
#'ppmdat <- ppmData(npoints = 10000,presences=presences, window = preds[[1]], covariates = preds)
#'species_formula <- presence ~ poly(X,2) + poly(Y,2) + poly(max_temp_hottest_month,2) + poly(annual_mean_precip,2) + poly(annual_mean_temp,2) + poly(distance_from_main_roads,2)
#'ft.ppm1 <- cvfit(species_formula = sp_form, object=ppmdat)
#'}
cvfit <- function(species_formula, bias_formula, object, method, control, nsim, ncores, ...){
  UseMethod("cvfit", object)
}

#' @export
cvfit.ppmData <- function(species_formula = presence/weights ~ 1,
                     bias_formula = NULL,
                     object,
                     method=c("lasso","glm","gam","ridge","ppmlasso"),
                     control=list(n.fit=20),
                     nsim=5,
                     ncores=1,...){

    ## set up a default method
    method <- match.arg(method)

    #need update the quadrature per-fit
    cvdatasets <- list()
    for(ii in seq_len(nsim)){
      cvdatasets[[ii]] <- pThin(ppmdata,p=p)
    }

    # cvdatasets <- lapply(1:nsim,pThin.ppmdat,ppmdata,p=p)
    cvppmfits <- list()
    for(ii in seq_len(nsim)){
      cvppmfits <- ppmFit(species_formula = species_formula,
                         bias_formula = bias_formula,
                         ppmdata = cvdatasets[[ii]]$train,
                         method = method)
    }

    return(list(cvppmfits,cvdatasets))

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
#' @export
pThin <- function(object, p = 0.25, seed=NULL, ...){
  UseMethod("cvfit", object)
}

#'@export
pThin.ppmData <- function(object, p = 0.25, seed=NULL){

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

   ppmdat_train$presences

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


#' @name blockSample
#' @description Do a block sample of sites
#' @param p The thinning probability (Bernoulli variable).
#' @author Skipton Woolley
#' @examples
#' library(ppmData)
#' presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
#' blockSample(presences[,-3],p.blocks=0.25,block.res=0.1)

blockSample <- function(sites, coords=c("X","Y"), p.blocks = 0.25, block.res = 1, seed = NULL){


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

  #locations of thinned points
  block_train=sites[-block.sites,]
  block_test=sites[block.sites,]

  # return as list
  res <- list(block_train_sites=block_train,
              block_test_sites=block_test)

  return(res)
}
