testthat::test_that('ppmFit lasso', {

  library(testthat)
  library(ppmData)
  library(ppmFit)
  library(terra)
  path <- system.file("extdata", package = "ppmData")
  lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
  preds <- rast(lst)
  presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
  npoints <- 1000

  # test with just pres, window & covars
  pd <- ppmData(presences=presences, window = preds[[1]], covariates = preds, npoints = npoints)

  # with species formula
  sp_form <- presence ~ poly(X,2) + poly(Y,2) + poly(max_temp_hottest_month,2) + poly(annual_mean_precip,2) + poly(annual_mean_temp,2) + poly(distance_from_main_roads,2)
  bias_form <- NULL
  ft.lasso1 <- ppmFit(species_formula = sp_form,bias_formula = bias_form, ppmdata=pd, method='lasso')
  expect_s3_class(ft.lasso1,"ppmFit")

  # with species + bias formula
  bias_form <- ~ -1 + distance_from_main_roads
  ft.lasso2 <- ppmFit(species_formula = sp_form, bias_formula = bias_form, ppmdata=pd, method='lasso')
  expect_s3_class(ft.lasso2,"ppmFit")

})
