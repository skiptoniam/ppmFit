## coefs
#'@importFrom glmnet coef.glmnet
coef.ppmFit <- function(object, s = NULL, exact = FALSE,...){
  coef.glmnet(object$ppm, s, exact, ...)
}
