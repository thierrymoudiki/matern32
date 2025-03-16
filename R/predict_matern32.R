#' Predict from Matern 3/2 model
#'
#' @param fit_obj Fitted object from \code{fit_matern32}
#' @param newx New data matrix
#' @param ci Not implemented yet
#'
#' @return
#' @export
#'
#' @examples
#' 
#' n <- 100 ; p <- 4
#' 
#' set.seed(456)
#' X <- matrix(rnorm(n * p), n, p) # no intercept!
#' y <- rnorm(n)
#' 
#' lams <- 10^seq(-10, 10, length.out = 50)
#' 
#' fit_obj <- fit_matern32(x = X, y = y, lambda = lams)
#' 
#' df <- data.frame(predict(fit_obj, X) - y)
#' colnames(df) <- paste0(round(lams, 2))
#' summary(df)
#' boxplot(df[, c(1, 10, 25, 35, 50)], 
#' main = "distribution of in sample bias", 
#' xlab = "lambda", ylab = "y_hat - y")
#' 
predict.matern32 <- function(fit_obj, newx, ci = NULL)
{
  if (is.vector(newx))
    newx <- t(newx)
  
  if (fit_obj$centering) 
  {
    # response was centered 
    if (!is.null(fit_obj$with_kmeans))
    {
      K_star <- matern32_kxstar_cpp(newx = as.matrix(matern32::my_scale(x = newx,
                                                                        xm = as.matrix(fit_obj$xm),
                                                                        xsd = as.vector(fit_obj$scales))), 
                                    x = fit_obj$scaled_x_clust, 
                                    l = fit_obj$l)
      
      return(drop(crossprod(K_star%*%fit_obj$coef)) + fit_obj$ym)
    } else {
      
      K_star <- matern32_kxstar_cpp(newx = as.matrix(matern32::my_scale(x = newx,
                                                                        xm = as.vector(fit_obj$xm),
                                                                        xsd = as.vector(fit_obj$scales))), 
                                    x = as.matrix(fit_obj$scaled_x), 
                                    l = as.vector(fit_obj$l))
      
      return(drop(crossprod(K_star, fit_obj$coef)) + fit_obj$ym)
    }
    
  } else {  
    
    # response wasn't centered 
    if (fit_obj$with_kmeans)
    {
      K_star <- matern32_kxstar_cpp(newx = as.matrix(matern32::my_scale(x = newx,
                                                                        xm = as.vector(fit_obj$xm),
                                                                        xsd = as.vector(fit_obj$scales))), 
                                    x = fit_obj$scaled_x_clust, 
                                    l = fit_obj$l)
      
      return(drop(crossprod(K_star%*%fit_obj$coef)))
    } else {
      
      K_star <- matern32_kxstar_cpp(newx = as.matrix(matern32::my_scale(x = newx,
                                                                        xm = as.vector(fit_obj$xm),
                                                                        xsd = as.vector(fit_obj$scales))), 
                                    x = fit_obj$scaled_x, 
                                    l = fit_obj$l)
      
      return(drop(K_star%*%fit_obj$coef))
    }
    
  }
  
}
#predict_matern32 <- memoise::memoize(predict_matern32)
