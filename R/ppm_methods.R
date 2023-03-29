## coefs
#'@importFrom glmnet coef.glmnet
coef.ppmFit <- function(object, s = NULL, exact = FALSE,...){
  coef.glmnet(object$ppm, s, exact, ...)
}


# loglikes
logLik.ppmFit <- function(object, ...){

  object$ppm$a0

  object$ppm$lambda

  mu <- predict(object)

  wts <- object$titbits$wts



  }


# if (family$family == "poisson")
# {
#   like = sum(ob.wt*(y*log(mu) - mu)) - sum(log(1:sum(y > 0))) - sum(alpha * as.vector(lambda)*abs(beta)) - sum(0.5 * (1 - alpha) * as.vector(lambda)*beta^2)
# }
# if (family$family == "binomial")
# {
#   like = sum(ob.wt*(y*log(mu) + (1 - y)*log(1 - mu))) - sum(as.vector(lambda)*abs(beta))
# }
