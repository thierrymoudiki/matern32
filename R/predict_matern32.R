#' Predict from Matern 3/2 model
#'
#' @param fit_obj 
#' @param newx 
#' @param ci 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' n <- 10 ; p <- 4
#' 
#' set.seed(456)
#' X <- matrix(rnorm(n * p), n, p) # no intercept!
#' y <- rnorm(n)
#' 
#' lams <- 10^seq(-5, 4, length.out = 50)
#' 
#' fit_obj <- fit_matern32(x = X, y = y, lambda = lams)
#' 
#' df <- data.frame(predict_matern32(fit_obj, newx = X) - y)
#' colnames(df) <- paste0(round(lams, 2))
#' summary(df)
#' boxplot(df[, c(1, 10, 25, 35, 50)], 
#' main = "distribution of bias", 
#' xlab = "lambda", ylab = "y_hat - y")
#' 
predict_matern32 <- function(fit_obj, newx, ci = NULL)
{
  if (is.vector(newx))
    newx <- t(newx)
  
  K_star <- matern32_kxstar_cpp(newx = as.matrix(matern32::my_scale(x = newx,
                                                                    xm = as.vector(fit_obj$xm),
                                                                    xsd = as.vector(fit_obj$scales))), 
                                x = fit_obj$scaled_x, 
                                l = fit_obj$l)
  
  return(drop(crossprod(K_star, fit_obj$coef)) + fit_obj$ym)
}
predict_matern32 <- memoise::memoize(predict_matern32)
