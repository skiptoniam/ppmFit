% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ppm_cv.R
\name{cvFit}
\alias{cvFit}
\title{A function for doing a very basic cross validation fitting based on a ppmData data object}
\usage{
cvFit(
  species_formula = presence/weights ~ 1,
  bias_formula = NULL,
  ppmdata,
  family = c("poisson", "binomial"),
  link = c("default", "log", "logit", "cloglog"),
  method = c("lasso", "ridge"),
  control = list(),
  titbits = TRUE,
  type = c("thin", "block"),
  nsim = 5,
  p = 0.25,
  res = 1,
  seed = NULL,
  ...
)
}
\arguments{
\item{species_formula}{The ppm model formula.}

\item{bias_formula}{Default is NULL. The idea will be to implement this as part of an integrated model.
Currently its just a nice way to keep track of which covariates are used in species or bias variables.
If you include this formula at the moment it will be merged into the species formula.}

\item{ppmdata}{ppmData object}

\item{family}{Character What family to use, current options are "poisson" for a poisson point process and "binomial" of an infinitely Weighted Logistic Regression (IWLR).}

\item{link}{Character What link function to use? Uses "default", "log" for poisson or "logit" for binomial. A user can specific "clogclog" for a binomial model.}

\item{method}{A method to fit the a ppm. Default is 'lasso', the alternative option is 'ridge' for ridge regression.}

\item{control}{list to pass to fitting functions.}

\item{titbits}{A boolean call to assess if you want extra model objects to be returned.}

\item{type}{character Use a "thin" approch or a "block" approach for setting up a cross-validation}

\item{nsim}{integer The number of cross validations to do (default is five).}

\item{p}{numeric The thinning probability (0,1).}

\item{res}{numeric The resolution of blocks if you want to try block resample}

\item{seed}{integer An integer to uses for set.seed(seed)}

\item{\\dots}{other things, not really used.}
}
\description{
This is a wrapper function which will fit a ppm model using a range of existing ppm fitting tools.
}
\examples{
\dontrun{
library(ppmData)
library(terra)
path <- system.file("extdata", package = "ppmData")
lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
preds <- rast(lst)
presences <- subset(snails,SpeciesID \%in\% "Tasmaphena sinclairi")
ppmdat <- ppmData(npoints = 10000,presences=presences, window = preds[[1]],
covariates = preds)
sp_form <- presence ~ poly(X,2) + poly(Y,2) +
poly(max_temp_hottest_month,2) + poly(annual_mean_precip,2) +
poly(annual_mean_temp,2) + poly(distance_from_main_roads,2)
cvft.ppm <- cvFit(species_formula = sp_form, ppmdata=ppmdat, family = "binomial", link = "logit")
}
}
\author{
Skipton Woolley
}
