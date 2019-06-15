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
  ym <- mean(y)
  scales <- my_scale(X)
  scaled_X <- scales$res
  xm <- scales$xm
  
  centered_y <- y - ym
  x1 <- scaled_X[,1]
  x2 <- scaled_X[,2]
  
  if (degree == 1)
  {
    XX <- cbind(scaled_X, x1*x2)
  }
  
  if (degree == 2)
  {
    XX <- cbind(scaled_X, scaled_X^2, 
          x1*x2)
  }
  
  if (degree == 3)
  {
    XX <- cbind(scaled_X, scaled_X^2, scaled_X^3, 
          x1*x2, x1*(x2^2), (x1^2)*x2, 
          (x1^2)*(x2^2))
  }
  
  return(list(fit_obj = list(coefficients = MASS::ginv(XX)%*%centered_y, 
                             degree = degree), 
              ym = ym, xm = xm))
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
  rescaled_X <- as.matrix(my_scale(newx, xm = object$xm))
  x1 <- rescaled_X[,1]
  x2 <- rescaled_X[,2]
  
  degree <- object$fit_obj$degree
  
  if (degree == 1)
  {
    XX <- cbind(rescaled_X, x1*x2)
  }
  
  if (degree == 2)
  {
    XX <- cbind(rescaled_X, rescaled_X^2, 
          x1*x2) 
  }
  
  if (degree == 3)
  {
    XX <- cbind(rescaled_X, rescaled_X^2, rescaled_X^3, 
          x1*x2, x1*(x2^2), (x1^2)*x2, 
          (x1^2)*(x2^2))
  }
  
  return(drop(object$ym + XX%*%as.numeric(object$fit_obj$coefficients)))
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
