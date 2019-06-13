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
fit_poly2 <- function(X, y)
{
  stopifnot(dim(X)[2] == 2)
  ym <- mean(y)
  scales <- my_scale(X)
  scaled_X <- scales$res
  xm <- scales$xm
  
  centered_y <- y - ym
  x1 <- scaled_X[,1]
  x2 <- scaled_X[,2]
  XX <- cbind(scaled_X, scaled_X^2, 
              x1*x2)
  
  return(list(fit_obj = .lm.fit(x = XX, y = centered_y), 
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
predict_poly2 <- function(object, newx)
{
  stopifnot(dim(newx)[2] == 2)
  rescaled_X <- as.matrix(my_scale(newx, xm = object$xm))
  x1 <- rescaled_X[,1]
  x2 <- rescaled_X[,2]
  XX <- cbind(rescaled_X, rescaled_X^2, 
              x1*x2)
  
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
