# Test poly2 --------------------------------------------------------------

#' Title
#'
#' @param X 
#' @param y 
#'
#' @return
#' @export
#'
#' @examples
fit_poly <- function(X, y, degree = 1)
{
  
  stopifnot(dim(X)[2] == 2)
  
  x1 <- X[, 1]
  x2 <- X[, 2]
  
  if (degree == 1)
  {
    XX <- cbind(X, x1*x2)
  }
  
  if (degree == 2)
  {
    XX <- cbind(X, X^2, x1*x2)
  }
  
  if (degree == 3)
  {
    XX <- cbind(X, X^2, X^3, 
          x1*x2, x1*(x2^2), (x1^2)*x2, 
          (x1^2)*(x2^2))
  }
  
  df <- data.frame(y = y, x = XX)
  
  obj <- MASS::lm.ridge(y ~ ., data=df, 
                        lambda = 10^seq(from=-10, to=10, 
                                        length.out = 100))

  obj$degree <- degree
  
  return(obj)
}


#' Title
#'
#' @param object 
#' @param newx 
#'
#' @return
#' @export
#'
#' @examples
predict_poly <- function(object, newx)
{
  stopifnot(dim(newx)[2] == 2)
  
  x1 <- newx[,1]
  x2 <- newx[,2]
  
  if (object$degree == 1)
  {
    XX <- cbind(X, x1*x2)
  }
  
  if (object$degree == 2)
  {
    XX <- cbind(X, X^2, x1*x2) 
  }
  
  if (object$degree == 3)
  {
    XX <- cbind(X, X^2, X^3, 
          x1*x2, x1*(x2^2), (x1^2)*x2, 
          (x1^2)*(x2^2))
  }
  
  scaled_x <- scale(XX, center = object$xm, scale = object$scales)
  
  return(drop(object$ym + scaled_x%*%as.numeric(object$coef[, which.min(object$GCV)])))
}


# n <- 25 ; p <- 2
# X <- matrix(rnorm(n * p), n, p) # no intercept!
# y <- rnorm(n)
# fit_obj <- fit_poly2(X, y)
# preds <- predict_poly2(fit_obj, X)
# 
# 
# plot(y, type = 'l', col = "blue")
# lines(preds, col = "red")
# 
