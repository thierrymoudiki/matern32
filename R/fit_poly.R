# Test poly2 --------------------------------------------------------------

#' Polynomial regression
#'
#' @param x Matrix of predictors
#' @param y Vector of response values 
#' @param degree Degree of the polynomial
#'
#' @return
#' @export
#'
#' @examples
fit_poly <- function(x, y, degree = 1)
{
  X <- as.matrix(x)
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
  
  class(obj) <- "poly"
  
  return(obj)
}


#' Predict method for polynomial regression
#'
#' @param object Fitted object
#' @param newx New data
#'
#' @return
#' @export
#'
#' @examples
predict.poly <- function(object, newx)
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
