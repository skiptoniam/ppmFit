# library(ppmData)
# library(terra)
# library(spatstat)
# library(glmnet)
# path <- system.file("extdata", package = "ppmData")
# lst <- list.files(path=path,pattern='*.tif',full.names = TRUE)
# preds <- rast(lst)
# presences <- subset(snails,SpeciesID %in% "Tasmaphena sinclairi")
#
# ppmdata <- ppmData(npoints = 100000,presences=presences,
#                   window = preds[[1]],
#                   covariates = preds)
#
# form <- ~ poly(max_temp_hottest_month, degree = 2) +
#   poly(annual_mean_precip, degree = 2) +
#   poly(annual_mean_temp, degree = 2) +
#   poly(distance_from_main_roads, degree = 2)
#
# # adaptive_lasso(form,ppmdata,"poisson")


adaptive_lasso <- function(form, ppmdata, family){

  ss <- as.spatstat(ppmdata)
  Q <- ss$Q
  ss.form <- update(form,Q~.)
  ss.fit <- ppm(ss.form,data=ss.covars,Poisson())

  ## set up the data structures for glmnet
  dat <- ss.fit$internal$glmfit$data
  f <- ss.fit$internal$glmfit$formula
  tms <- ss.fit$internal$glmfit$terms
  yy <- dat[,2]
  n <- sum(ppmdata$ppmData$presence)
  N <-  length(yy)
  wts <- dat[,1]
  mf <- model.frame(formula = f, data = dat) # weights need to be name of weights in ppp
  mt <- terms(mf)
  Q.new <- model.matrix(mt,mf)

  ## need to calculate the sigma and hessian of the inital ppm.fit
  H <- vcov(ss.fit, what="fisher", hessian=TRUE, fine=FALSE)
  Sigma <- vcov(ss.fit,fine=FALSE)

  ## Initial lasso - really just to set up the lambdas
  fit.lasso <- glmnet(Q.new[,-1], yy, family="poisson", alpha=1, weights=wts,
                      penalty.factor=c(rep(sum(wts),number.covariates-1)))
  coef.theta.lasso <- as.matrix(rbind(fit.lasso$a0, fit.lasso$beta))
  xi.lasso <- as.matrix(fit.lasso$lambda)

  ## Initial estimate
  num.covars <- nrow(H)
  tmp <- c(rep(0,num.covars-1)) # drop intercept
  theta.init <- tmp+runif(length(tmp),-.1,.1)

  ## Adaptive Lasso
  fit.alasso <- glmnet(Q.new[,-1], yy, family="poisson", alpha=1, weights=wts,
                       penalty.factor=c(sum(wts)/abs(theta.init)))
  coef.theta.alasso <- as.matrix(rbind(fit.alasso$a0,fit.alasso$beta))
  xi.alasso <- as.matrix(fit.alasso$lambda)


  optim.inf.criteria <- 1e30
  index.optim.theta <- NULL
  inf.criteria <- NULL
  for (i in 1:length(xi.lasso)) {
    theta <- coef.theta.lasso[,i]
    ind.zero <- which(theta==0)
    logpseudolike <- sum(wts*(yy*(Q.new%*%theta)-exp(Q.new%*%theta)))
    Sigma.new <- Sigma
    H.new <- H
    if (length(ind.zero)!=0) {
      Sigma.new <- Sigma.new[-ind.zero,-ind.zero]
      H.new <- H.new[-ind.zero,-ind.zero]
    }
    ## cBIC criterion
    inf.criteria[i] <- -2*logpseudolike + (sum(diag(H.new%*%Sigma.new)))*log(n)
    if (inf.criteria[i] < optim.inf.criteria) {
      index.optim.theta <- i
      optim.inf.criteria <- inf.criteria[i]
    }
  }
  theta.lasso.cBIC <- coef.theta.lasso[,index.optim.theta]
  xi.opt.lasso.cBIC <- xi.lasso[index.optim.theta]
  xi.opt.lasso.cBIC # optimal xi parameter

  optim.inf.criteria <- 1e30
  index.optim.theta <- NULL
  inf.criteria <- NULL

  for (i in 1:length(xi.alasso)) {
    theta <- coef.theta.alasso[,i]
    ind.zero <- which(theta==0)
    logpseudolike <- sum(wts*(yy*(Q.new%*%theta)-exp(Q.new%*%theta)))
    Sigma.new <- Sigma
    H.new <- H
    if (length(ind.zero)!=0) {
      Sigma.new <- Sigma.new[-ind.zero,-ind.zero]
      H.new <- H.new[-ind.zero,-ind.zero]
    }
    ## cERIC criterion
    inf.criteria[i] <- -2*logpseudolike + (sum(diag(H.new%*%Sigma.new)))*log(n/(N*xi.alasso[i]))
    if (inf.criteria[i] < optim.inf.criteria) {
      index.optim.theta <- i
      optim.inf.criteria <- inf.criteria[i]
    }
  }
  theta.alasso.cERIC <- coef.theta.alasso[,index.optim.theta]
  xi.opt.alasso.cERIC <- xi.alasso[index.optim.theta]
  xi.opt.alasso.cERIC # optimal xi parameter


}



