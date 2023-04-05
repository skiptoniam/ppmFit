## coefs
#'@importFrom glmnet coef.glmnet
coef.ppmFit <- function(object, s = NULL, exact = FALSE,...){
  coef.glmnet(object$ppm, s, exact, ...)
}


# loglikes
logLik.ppmFit <- function(object, ...){

  object$ppm$a0

  object$ppm$lambda

  object$ppm$beta

  mu <- as.numeric(predict(object),quad.only=FALSE)
  y <- as.numeric(object$titbits$y)
  wts <- as.numeric(object$titbits$wts)

  if (family$family == "poisson"){
      like = sum(wts*(y*log(mu) - mu)) - sum(log(1:sum(y > 0))) - sum(alpha * as.vector(lambda)*abs(beta)) - sum(0.5 * (1 - alpha) * as.vector(lambda)*beta^2)
  }
  if (family$family == "binomial"){
      like = sum(wts*(y*log(mu) + (1 - y)*log(1 - mu)),na.rm = TRUE) - sum(as.vector(lambda)*abs(beta))
  }

}
