testthat::test_that('ppmFit lasso', {

  library(testthat)
  library(ppmData)
  library(ppmFit)
  path <- system.file("extdata", package = "ppmData")
  lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
  preds <- rast(lst)
  presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
  npoints <- 1000

  # test with just pres, window & covars
  pd <- ppmData(presences=presences, window = preds[[1]], covariates = preds, npoints = npoints)

  # with species formula
  sp_form <- presence ~  poly(max_temp_hottest_month,2) + poly(annual_mean_precip,2) + poly(annual_mean_temp,2) + poly(distance_from_main_roads,2)
  bias_form <- NULL
  ft.lasso1 <- ppmFit(species_formula = sp_form,bias_formula = bias_form, ppmdata=pd)

  ## test self predict at observer points
  pred <- predict(ft.lasso1)
  expect_type(pred,"double")

  ## test predict to just quadrature sites
  pred.quad <- predict(ft.lasso1,quad.only=TRUE)
  expect_type(pred.quad,"double")

  ## test predict to SpatRaster - terra
  pred.rast <- predict(ft.lasso1, newdata=preds,quad.only=TRUE)
  expect_s4_class(pred.rast,"SpatRaster")


})
